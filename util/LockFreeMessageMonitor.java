package drake.util;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import lcm.lcm.LCM;
import lcm.lcm.LCMDataInputStream;
import lcm.lcm.LCMSubscriber;

public class LockFreeMessageMonitor implements LCMSubscriber {

  private final Map<String, ConcurrentCopier<LCMMessageData>> copiers = new ConcurrentHashMap<String, ConcurrentCopier<LCMMessageData>>();
  private final Map<String, LCMMessageData> last_returned_data = new HashMap<String, LCMMessageData>();
  
  public LockFreeMessageMonitor() {
  }

  @Override
  public void messageReceived(LCM lcm, String channel, LCMDataInputStream ins) {
    ConcurrentCopier<LCMMessageData> copier = copiers.get(channel);
    if (copier == null) {
      copier = new ConcurrentCopier<LockFreeMessageMonitor.LCMMessageData>(new LCMMessageDataBuilder());
      copiers.put(channel, copier);
    }
    
    LCMMessageData data = copier.getCopyForWriting();
    int available = ins.available();
    if (data.byte_array == null || data.byte_array.length != available) { // TODO: is it OK to have a byte array that is too large?
      data.byte_array = new byte[available];
    }
    try {
      ins.readFully(data.byte_array);
    } catch (IOException e) {
      System.err.println("MultipleMessageMonitor exception on channel " + channel);
      e.printStackTrace();
    }
    copier.commit();
  }

  /*
   * Assumed to be called from a single thread
   */
  public Map<String, byte[]> getMessages() {
    Map<String, byte[]> ret = new HashMap<String, byte[]>();
    for (String channel : copiers.keySet()) {
      ConcurrentCopier<LCMMessageData> copier = copiers.get(channel);
      LCMMessageData data = copier.getCopyForReading();
      if (data != null && data != last_returned_data.get(channel)) {
        ret.put(channel, Arrays.copyOf(data.byte_array, data.byte_array.length));
        last_returned_data.put(channel, data);
      }
    }
    
    return ret;
  }

  private static class LCMMessageData {
    public byte[] byte_array;
  }

  private static class LCMMessageDataBuilder implements Builder<LCMMessageData> {

    @Override
    public LCMMessageData newInstance() {
      return new LCMMessageData();
    }
  }
}

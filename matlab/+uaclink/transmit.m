function [waveform, info] = transmit(payload, sequence, frameType, cfg)
%TRANSMIT Serialize, CRC, and modulate a wire frame including HFM and pilots.
[bytes,bits]=uaclink.encode_frame(payload,sequence,frameType,cfg);
[waveform,layout]=uac2.make_packet(bits,cfg);
info=struct('wire_bytes',bytes,'n_bits',numel(bits), ...
    'duration_s',layout.packet_duration_s,'sequence',sequence,'type',frameType);
end

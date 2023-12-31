OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(2.2709742) q[0];
rz(2.119996) q[1];
sx q[1];
rz(-2.8586913) q[1];
sx q[1];
rz(0.14970782) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5538841) q[0];
sx q[0];
rz(-2.2061078) q[0];
sx q[0];
rz(-2.5250838) q[0];
rz(2.6354191) q[2];
sx q[2];
rz(-1.0867599) q[2];
sx q[2];
rz(-1.5278221) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40587273) q[1];
sx q[1];
rz(-2.957203) q[1];
sx q[1];
rz(-1.7806446) q[1];
rz(1.2655067) q[3];
sx q[3];
rz(-0.57170924) q[3];
sx q[3];
rz(-1.3523462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(0.70409888) q[2];
rz(0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(0.63252226) q[0];
rz(2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(0.65223637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1855542) q[0];
sx q[0];
rz(-1.5421876) q[0];
sx q[0];
rz(-1.5584598) q[0];
x q[1];
rz(1.2781497) q[2];
sx q[2];
rz(-1.3424982) q[2];
sx q[2];
rz(1.6739068) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10416874) q[1];
sx q[1];
rz(-1.2836604) q[1];
sx q[1];
rz(-1.3377405) q[1];
x q[2];
rz(-0.60450508) q[3];
sx q[3];
rz(-1.0009871) q[3];
sx q[3];
rz(-1.9250211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(2.2913349) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(1.7920378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3916546) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(-1.2664938) q[0];
rz(-pi) q[1];
rz(1.6367958) q[2];
sx q[2];
rz(-1.5335576) q[2];
sx q[2];
rz(0.50022349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.459356) q[1];
sx q[1];
rz(-0.56486928) q[1];
sx q[1];
rz(0.45046803) q[1];
x q[2];
rz(-1.5924256) q[3];
sx q[3];
rz(-0.4261741) q[3];
sx q[3];
rz(0.91859761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0321908) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(0.90399495) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.7310463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(0.036380336) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0050126652) q[0];
sx q[0];
rz(-2.4653325) q[0];
sx q[0];
rz(-1.9146634) q[0];
rz(1.245001) q[2];
sx q[2];
rz(-2.1278283) q[2];
sx q[2];
rz(-1.4404802) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.034678) q[1];
sx q[1];
rz(-2.2639096) q[1];
sx q[1];
rz(0.17130674) q[1];
x q[2];
rz(-0.23457228) q[3];
sx q[3];
rz(-1.4751225) q[3];
sx q[3];
rz(-0.58805874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2146384) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(0.67048091) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(2.9751076) q[0];
rz(-2.3855551) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(-2.9072445) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4011824) q[0];
sx q[0];
rz(-0.64490841) q[0];
sx q[0];
rz(-2.5572204) q[0];
rz(-3.0460065) q[2];
sx q[2];
rz(-2.0628953) q[2];
sx q[2];
rz(1.0620067) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6926596) q[1];
sx q[1];
rz(-1.7198791) q[1];
sx q[1];
rz(-2.3922608) q[1];
rz(-1.3393199) q[3];
sx q[3];
rz(-0.33208648) q[3];
sx q[3];
rz(1.0486697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70871893) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(0.22053545) q[2];
rz(2.7045414) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(-2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-0.81714001) q[0];
rz(0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(1.9979427) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8457984) q[0];
sx q[0];
rz(-2.7685071) q[0];
sx q[0];
rz(-0.11008115) q[0];
x q[1];
rz(2.8191889) q[2];
sx q[2];
rz(-0.69903261) q[2];
sx q[2];
rz(-2.2088745) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3559349) q[1];
sx q[1];
rz(-0.4520843) q[1];
sx q[1];
rz(-1.362734) q[1];
rz(-pi) q[2];
rz(-1.2915217) q[3];
sx q[3];
rz(-2.0719299) q[3];
sx q[3];
rz(-0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5722826) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(1.8547159) q[0];
rz(-0.6634179) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(1.2333262) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7147303) q[0];
sx q[0];
rz(-0.38422248) q[0];
sx q[0];
rz(-1.5412488) q[0];
rz(-2.2912824) q[2];
sx q[2];
rz(-1.0734953) q[2];
sx q[2];
rz(0.26956272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8804633) q[1];
sx q[1];
rz(-1.556734) q[1];
sx q[1];
rz(-0.1111828) q[1];
rz(0.4864278) q[3];
sx q[3];
rz(-1.8842116) q[3];
sx q[3];
rz(0.13669554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.037501637) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0163517) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(2.9267172) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.8833556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5364752) q[0];
sx q[0];
rz(-2.0015171) q[0];
sx q[0];
rz(-2.9933948) q[0];
rz(-2.3636742) q[2];
sx q[2];
rz(-0.81233835) q[2];
sx q[2];
rz(1.1493491) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59626034) q[1];
sx q[1];
rz(-0.55570554) q[1];
sx q[1];
rz(0.19821367) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85556742) q[3];
sx q[3];
rz(-1.1547525) q[3];
sx q[3];
rz(0.95203979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68226472) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11583081) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(-2.8905706) q[0];
rz(-2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(-0.02773157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6513034) q[0];
sx q[0];
rz(-0.9285183) q[0];
sx q[0];
rz(-2.519033) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8856144) q[2];
sx q[2];
rz(-2.417649) q[2];
sx q[2];
rz(-2.1997423) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30222826) q[1];
sx q[1];
rz(-1.4711958) q[1];
sx q[1];
rz(-2.539413) q[1];
x q[2];
rz(-0.025891993) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(0.47749146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1256844) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(1.998418) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(-2.4699396) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-0.25751105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6637708) q[0];
sx q[0];
rz(-1.1241233) q[0];
sx q[0];
rz(-1.1429943) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72980482) q[2];
sx q[2];
rz(-0.68241718) q[2];
sx q[2];
rz(-2.6953816) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.820206) q[1];
sx q[1];
rz(-0.69637978) q[1];
sx q[1];
rz(-0.022298261) q[1];
rz(-2.8994843) q[3];
sx q[3];
rz(-1.1032411) q[3];
sx q[3];
rz(-0.20413354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8132849) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(-0.47484067) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3474779) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-2.4702934) q[2];
sx q[2];
rz(-0.57325267) q[2];
sx q[2];
rz(-0.43043955) q[2];
rz(2.0444617) q[3];
sx q[3];
rz(-0.53285014) q[3];
sx q[3];
rz(0.22900029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

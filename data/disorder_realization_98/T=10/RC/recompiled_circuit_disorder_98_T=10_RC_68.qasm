OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(-2.1455278) q[0];
sx q[0];
rz(-2.2709742) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(-0.14970782) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4275442) q[0];
sx q[0];
rz(-0.85435003) q[0];
sx q[0];
rz(2.2358405) q[0];
x q[1];
rz(2.6354191) q[2];
sx q[2];
rz(-2.0548327) q[2];
sx q[2];
rz(-1.6137705) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61923164) q[1];
sx q[1];
rz(-1.751096) q[1];
sx q[1];
rz(3.1027604) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2655067) q[3];
sx q[3];
rz(-0.57170924) q[3];
sx q[3];
rz(-1.7892464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8831138) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(0.70409888) q[2];
rz(2.1885833) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(2.5090704) q[0];
rz(2.6951492) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(-0.65223637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61440496) q[0];
sx q[0];
rz(-1.5584649) q[0];
sx q[0];
rz(3.1129818) q[0];
rz(-pi) q[1];
rz(-0.23806439) q[2];
sx q[2];
rz(-1.2859584) q[2];
sx q[2];
rz(-0.1711947) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7420885) q[1];
sx q[1];
rz(-1.3474476) q[1];
sx q[1];
rz(0.29466596) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2324123) q[3];
sx q[3];
rz(-2.069807) q[3];
sx q[3];
rz(3.1391075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(-1.1616421) q[2];
rz(1.15796) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(-2.2913349) q[0];
rz(2.6440874) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(-1.3495548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36950668) q[0];
sx q[0];
rz(-1.8096576) q[0];
sx q[0];
rz(2.4580965) q[0];
rz(-pi) q[1];
rz(3.1042728) q[2];
sx q[2];
rz(-1.5048426) q[2];
sx q[2];
rz(2.0685591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68223665) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(-0.45046803) q[1];
rz(-0.0098185929) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(-2.1992418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(-0.90399495) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(0.63181216) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(-3.1052123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.13658) q[0];
sx q[0];
rz(-2.4653325) q[0];
sx q[0];
rz(-1.2269292) q[0];
rz(-pi) q[1];
rz(-2.5600299) q[2];
sx q[2];
rz(-1.8459324) q[2];
sx q[2];
rz(3.0951701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2991043) q[1];
sx q[1];
rz(-0.71055382) q[1];
sx q[1];
rz(1.7732265) q[1];
x q[2];
rz(1.6691469) q[3];
sx q[3];
rz(-1.8042759) q[3];
sx q[3];
rz(-1.0055621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.1882163) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(2.9072445) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74041022) q[0];
sx q[0];
rz(-2.4966842) q[0];
sx q[0];
rz(0.58437225) q[0];
rz(-pi) q[1];
rz(1.0767897) q[2];
sx q[2];
rz(-1.4865808) q[2];
sx q[2];
rz(-2.6780724) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8824132) q[1];
sx q[1];
rz(-0.83173527) q[1];
sx q[1];
rz(1.773136) q[1];
x q[2];
rz(-1.2469532) q[3];
sx q[3];
rz(-1.4959335) q[3];
sx q[3];
rz(2.4002241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(-0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-0.81714001) q[0];
rz(-0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957942) q[0];
sx q[0];
rz(-2.7685071) q[0];
sx q[0];
rz(0.11008115) q[0];
rz(-pi) q[1];
rz(1.3104865) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(-1.7973763) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.027187849) q[1];
sx q[1];
rz(-1.6611551) q[1];
sx q[1];
rz(1.1272217) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51795824) q[3];
sx q[3];
rz(-1.3266139) q[3];
sx q[3];
rz(-1.9675919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(-2.9135381) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(2.4781748) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(-1.9082665) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42686233) q[0];
sx q[0];
rz(-2.7573702) q[0];
sx q[0];
rz(1.5412488) q[0];
rz(-pi) q[1];
rz(0.8823422) q[2];
sx q[2];
rz(-2.2922278) q[2];
sx q[2];
rz(1.7989858) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8804633) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(-0.1111828) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60621467) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(-1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.037501637) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(-1.0160149) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(1.4455147) q[0];
rz(-0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(-1.8833556) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90342605) q[0];
sx q[0];
rz(-1.4362207) q[0];
sx q[0];
rz(-2.0057136) q[0];
x q[1];
rz(2.2082981) q[2];
sx q[2];
rz(-1.0273232) q[2];
sx q[2];
rz(2.9535434) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3360577) q[1];
sx q[1];
rz(-1.6748669) q[1];
sx q[1];
rz(0.54688262) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6120841) q[3];
sx q[3];
rz(-2.2141075) q[3];
sx q[3];
rz(0.95638004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(2.5788467) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(3.1138611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3656317) q[0];
sx q[0];
rz(-2.2790475) q[0];
sx q[0];
rz(-2.2328949) q[0];
rz(-pi) q[1];
rz(0.8717732) q[2];
sx q[2];
rz(-1.36424) q[2];
sx q[2];
rz(-2.2733462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7293538) q[1];
sx q[1];
rz(-2.5322399) q[1];
sx q[1];
rz(-0.17462294) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.025891993) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(-2.6641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1256844) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(-2.9987191) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(2.8840816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2902381) q[0];
sx q[0];
rz(-2.5332753) q[0];
sx q[0];
rz(0.71382199) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4117878) q[2];
sx q[2];
rz(-2.4591755) q[2];
sx q[2];
rz(-2.6953816) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7911437) q[1];
sx q[1];
rz(-0.87462438) q[1];
sx q[1];
rz(1.5894366) q[1];
rz(-pi) q[2];
rz(-2.0503644) q[3];
sx q[3];
rz(-1.3551095) q[3];
sx q[3];
rz(1.2558162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(-2.5349687) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(2.9150302) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(-2.6735641) q[2];
sx q[2];
rz(-1.2266908) q[2];
sx q[2];
rz(-2.5897349) q[2];
rz(2.0541035) q[3];
sx q[3];
rz(-1.3369505) q[3];
sx q[3];
rz(-1.7575775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
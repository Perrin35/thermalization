OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.8340303) q[0];
sx q[0];
rz(-2.7355255) q[0];
sx q[0];
rz(1.9598444) q[0];
rz(1.740068) q[1];
sx q[1];
rz(-0.47456804) q[1];
sx q[1];
rz(0.13303703) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2849738) q[0];
sx q[0];
rz(-0.70278888) q[0];
sx q[0];
rz(-1.982559) q[0];
rz(0.51271397) q[2];
sx q[2];
rz(-1.1774592) q[2];
sx q[2];
rz(-0.31189298) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2084004) q[1];
sx q[1];
rz(-2.7504535) q[1];
sx q[1];
rz(-2.1927298) q[1];
rz(-2.1054793) q[3];
sx q[3];
rz(-2.4991192) q[3];
sx q[3];
rz(-1.5657263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1349606) q[2];
sx q[2];
rz(-1.0198318) q[2];
sx q[2];
rz(2.4572065) q[2];
rz(-2.0002174) q[3];
sx q[3];
rz(-2.8751774) q[3];
sx q[3];
rz(-0.086708955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619277) q[0];
sx q[0];
rz(-2.4517038) q[0];
sx q[0];
rz(2.6806114) q[0];
rz(-2.0224679) q[1];
sx q[1];
rz(-0.42367595) q[1];
sx q[1];
rz(0.31203312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082816007) q[0];
sx q[0];
rz(-3.1110373) q[0];
sx q[0];
rz(-2.729134) q[0];
rz(-pi) q[1];
x q[1];
rz(0.076032555) q[2];
sx q[2];
rz(-0.8113779) q[2];
sx q[2];
rz(-1.8234314) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65822433) q[1];
sx q[1];
rz(-1.7464356) q[1];
sx q[1];
rz(-3.1356922) q[1];
rz(-pi) q[2];
rz(2.8617763) q[3];
sx q[3];
rz(-0.93799611) q[3];
sx q[3];
rz(1.0029083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95388874) q[2];
sx q[2];
rz(-2.0362594) q[2];
sx q[2];
rz(2.7375431) q[2];
rz(-2.7814501) q[3];
sx q[3];
rz(-1.6481684) q[3];
sx q[3];
rz(-0.68296105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78524154) q[0];
sx q[0];
rz(-2.1389565) q[0];
sx q[0];
rz(-2.8175765) q[0];
rz(2.0815381) q[1];
sx q[1];
rz(-2.8476604) q[1];
sx q[1];
rz(2.4411328) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4361079) q[0];
sx q[0];
rz(-1.4027081) q[0];
sx q[0];
rz(-1.4072493) q[0];
rz(-0.21493463) q[2];
sx q[2];
rz(-1.2329458) q[2];
sx q[2];
rz(-0.95906729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1478657) q[1];
sx q[1];
rz(-0.7997457) q[1];
sx q[1];
rz(-0.32983671) q[1];
rz(2.9459314) q[3];
sx q[3];
rz(-1.1563468) q[3];
sx q[3];
rz(0.54563845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6274274) q[2];
sx q[2];
rz(-2.8943987) q[2];
sx q[2];
rz(-2.3449281) q[2];
rz(-1.944815) q[3];
sx q[3];
rz(-1.6007042) q[3];
sx q[3];
rz(0.59320199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45348039) q[0];
sx q[0];
rz(-2.1402833) q[0];
sx q[0];
rz(1.8185115) q[0];
rz(0.46609136) q[1];
sx q[1];
rz(-0.90704647) q[1];
sx q[1];
rz(-1.1493433) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2043484) q[0];
sx q[0];
rz(-2.4597628) q[0];
sx q[0];
rz(-2.3005062) q[0];
x q[1];
rz(2.7485899) q[2];
sx q[2];
rz(-2.1519024) q[2];
sx q[2];
rz(-0.71497922) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9117078) q[1];
sx q[1];
rz(-0.20710342) q[1];
sx q[1];
rz(-2.9207499) q[1];
x q[2];
rz(-1.2364616) q[3];
sx q[3];
rz(-1.4282246) q[3];
sx q[3];
rz(0.73745525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.484802) q[2];
sx q[2];
rz(-0.66978407) q[2];
sx q[2];
rz(2.316324) q[2];
rz(1.8615865) q[3];
sx q[3];
rz(-1.5214336) q[3];
sx q[3];
rz(0.72117225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44329062) q[0];
sx q[0];
rz(-1.8809141) q[0];
sx q[0];
rz(1.7753506) q[0];
rz(2.0514936) q[1];
sx q[1];
rz(-1.5049728) q[1];
sx q[1];
rz(-1.8025788) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96025673) q[0];
sx q[0];
rz(-0.68013817) q[0];
sx q[0];
rz(2.8614317) q[0];
rz(-pi) q[1];
rz(2.5191804) q[2];
sx q[2];
rz(-0.96660766) q[2];
sx q[2];
rz(2.6780918) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.73698649) q[1];
sx q[1];
rz(-1.6885719) q[1];
sx q[1];
rz(0.21356689) q[1];
x q[2];
rz(-0.33356922) q[3];
sx q[3];
rz(-1.8954479) q[3];
sx q[3];
rz(-2.1683104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1947386) q[2];
sx q[2];
rz(-1.7297435) q[2];
sx q[2];
rz(-2.8080688) q[2];
rz(-2.5528095) q[3];
sx q[3];
rz(-2.6684561) q[3];
sx q[3];
rz(0.80140448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81353417) q[0];
sx q[0];
rz(-1.7840339) q[0];
sx q[0];
rz(-1.300977) q[0];
rz(-1.0282372) q[1];
sx q[1];
rz(-2.6515549) q[1];
sx q[1];
rz(-0.44233826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16602941) q[0];
sx q[0];
rz(-1.1021758) q[0];
sx q[0];
rz(-1.8450849) q[0];
rz(2.6643941) q[2];
sx q[2];
rz(-2.0522567) q[2];
sx q[2];
rz(1.9419326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7117911) q[1];
sx q[1];
rz(-0.45908976) q[1];
sx q[1];
rz(0.38021047) q[1];
rz(-2.836727) q[3];
sx q[3];
rz(-0.33434048) q[3];
sx q[3];
rz(1.5078733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3698795) q[2];
sx q[2];
rz(-1.5481202) q[2];
sx q[2];
rz(0.4717353) q[2];
rz(-1.6996023) q[3];
sx q[3];
rz(-2.580018) q[3];
sx q[3];
rz(-0.26427463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1137375) q[0];
sx q[0];
rz(-1.6738482) q[0];
sx q[0];
rz(2.9748528) q[0];
rz(-0.79404229) q[1];
sx q[1];
rz(-1.200095) q[1];
sx q[1];
rz(0.099743191) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3152811) q[0];
sx q[0];
rz(-2.4874685) q[0];
sx q[0];
rz(-2.099206) q[0];
x q[1];
rz(1.0641618) q[2];
sx q[2];
rz(-0.68175137) q[2];
sx q[2];
rz(-3.0118941) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9236026) q[1];
sx q[1];
rz(-1.0031219) q[1];
sx q[1];
rz(-1.7204549) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2325956) q[3];
sx q[3];
rz(-1.2332877) q[3];
sx q[3];
rz(2.8765529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2793067) q[2];
sx q[2];
rz(-2.0760832) q[2];
sx q[2];
rz(0.12969895) q[2];
rz(1.8779523) q[3];
sx q[3];
rz(-1.2010776) q[3];
sx q[3];
rz(0.14550801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97067863) q[0];
sx q[0];
rz(-0.92249528) q[0];
sx q[0];
rz(-2.7150735) q[0];
rz(2.049394) q[1];
sx q[1];
rz(-2.009232) q[1];
sx q[1];
rz(-1.0027286) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9364081) q[0];
sx q[0];
rz(-2.4800088) q[0];
sx q[0];
rz(-1.9993881) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10522225) q[2];
sx q[2];
rz(-2.4159523) q[2];
sx q[2];
rz(-1.9228455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59834023) q[1];
sx q[1];
rz(-1.9045254) q[1];
sx q[1];
rz(-2.0327492) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9055481) q[3];
sx q[3];
rz(-1.1268643) q[3];
sx q[3];
rz(0.027396552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2670474) q[2];
sx q[2];
rz(-2.5858877) q[2];
sx q[2];
rz(3.0181273) q[2];
rz(-0.63001436) q[3];
sx q[3];
rz(-1.29653) q[3];
sx q[3];
rz(0.71648487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89710871) q[0];
sx q[0];
rz(-0.69863313) q[0];
sx q[0];
rz(-0.83936349) q[0];
rz(-0.15946236) q[1];
sx q[1];
rz(-2.1493252) q[1];
sx q[1];
rz(2.2119567) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63422257) q[0];
sx q[0];
rz(-2.3593724) q[0];
sx q[0];
rz(-0.36391516) q[0];
rz(-0.18126296) q[2];
sx q[2];
rz(-2.5428704) q[2];
sx q[2];
rz(-1.8229654) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8690455) q[1];
sx q[1];
rz(-2.2490977) q[1];
sx q[1];
rz(2.9806304) q[1];
x q[2];
rz(1.2655018) q[3];
sx q[3];
rz(-2.4809894) q[3];
sx q[3];
rz(-1.5098287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.9440426) q[2];
sx q[2];
rz(-1.9746747) q[2];
sx q[2];
rz(-0.76773947) q[2];
rz(-0.36457101) q[3];
sx q[3];
rz(-1.7578099) q[3];
sx q[3];
rz(-3.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1345054) q[0];
sx q[0];
rz(-0.61408478) q[0];
sx q[0];
rz(2.2579204) q[0];
rz(-0.20835749) q[1];
sx q[1];
rz(-0.50269428) q[1];
sx q[1];
rz(-1.3349894) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55157298) q[0];
sx q[0];
rz(-0.087170426) q[0];
sx q[0];
rz(0.7256736) q[0];
rz(2.8284188) q[2];
sx q[2];
rz(-2.5921481) q[2];
sx q[2];
rz(0.012481364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.523004) q[1];
sx q[1];
rz(-2.4135313) q[1];
sx q[1];
rz(0.87597202) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6649762) q[3];
sx q[3];
rz(-0.94175613) q[3];
sx q[3];
rz(-1.4339989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9779382) q[2];
sx q[2];
rz(-1.6341219) q[2];
sx q[2];
rz(0.046023544) q[2];
rz(2.9764791) q[3];
sx q[3];
rz(-2.8486227) q[3];
sx q[3];
rz(1.9273812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56275) q[0];
sx q[0];
rz(-3.0414707) q[0];
sx q[0];
rz(1.1895251) q[0];
rz(2.9764755) q[1];
sx q[1];
rz(-1.4585635) q[1];
sx q[1];
rz(-0.30452902) q[1];
rz(-0.44331917) q[2];
sx q[2];
rz(-1.8675315) q[2];
sx q[2];
rz(0.055589614) q[2];
rz(2.6611609) q[3];
sx q[3];
rz(-1.1050925) q[3];
sx q[3];
rz(2.1375124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

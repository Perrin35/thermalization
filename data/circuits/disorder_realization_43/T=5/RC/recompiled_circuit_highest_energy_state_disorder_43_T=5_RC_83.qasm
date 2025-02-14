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
rz(2.9146258) q[0];
sx q[0];
rz(-0.6114971) q[0];
sx q[0];
rz(-2.5831232) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(-1.7668626) q[1];
sx q[1];
rz(-3.0629646) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0315379) q[0];
sx q[0];
rz(-2.1221836) q[0];
sx q[0];
rz(-0.64161513) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4628381) q[2];
sx q[2];
rz(-1.4732483) q[2];
sx q[2];
rz(2.7979948) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3529571) q[1];
sx q[1];
rz(-0.26600263) q[1];
sx q[1];
rz(-2.8255844) q[1];
rz(-pi) q[2];
rz(1.2492248) q[3];
sx q[3];
rz(-1.235629) q[3];
sx q[3];
rz(-0.17091076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2416396) q[2];
sx q[2];
rz(-0.048154801) q[2];
sx q[2];
rz(-2.107035) q[2];
rz(-0.11040802) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(-0.30606562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13650525) q[0];
sx q[0];
rz(-1.445048) q[0];
sx q[0];
rz(-1.9553631) q[0];
rz(1.8738481) q[1];
sx q[1];
rz(-2.6823951) q[1];
sx q[1];
rz(1.3892106) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88436172) q[0];
sx q[0];
rz(-2.5116034) q[0];
sx q[0];
rz(1.0347741) q[0];
rz(-pi) q[1];
rz(-0.21135881) q[2];
sx q[2];
rz(-2.6508001) q[2];
sx q[2];
rz(0.15891128) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0276549) q[1];
sx q[1];
rz(-1.9009556) q[1];
sx q[1];
rz(0.3056674) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1950483) q[3];
sx q[3];
rz(-1.6485751) q[3];
sx q[3];
rz(0.7796703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17239751) q[2];
sx q[2];
rz(-2.0266504) q[2];
sx q[2];
rz(-1.4364852) q[2];
rz(1.2596333) q[3];
sx q[3];
rz(-0.45537046) q[3];
sx q[3];
rz(-3.1148541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7773975) q[0];
sx q[0];
rz(-1.7561678) q[0];
sx q[0];
rz(1.2170323) q[0];
rz(0.2150391) q[1];
sx q[1];
rz(-1.8937078) q[1];
sx q[1];
rz(-1.7020114) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0501409) q[0];
sx q[0];
rz(-0.9147075) q[0];
sx q[0];
rz(-1.4138282) q[0];
x q[1];
rz(-2.0820175) q[2];
sx q[2];
rz(-1.6523835) q[2];
sx q[2];
rz(-0.64543426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90104687) q[1];
sx q[1];
rz(-2.8786297) q[1];
sx q[1];
rz(2.0790624) q[1];
rz(-2.923449) q[3];
sx q[3];
rz(-1.8379595) q[3];
sx q[3];
rz(-1.657079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3415459) q[2];
sx q[2];
rz(-1.505082) q[2];
sx q[2];
rz(0.81614196) q[2];
rz(3.1331565) q[3];
sx q[3];
rz(-2.4237027) q[3];
sx q[3];
rz(1.5754023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4207882) q[0];
sx q[0];
rz(-1.9705462) q[0];
sx q[0];
rz(0.26914832) q[0];
rz(-1.3828269) q[1];
sx q[1];
rz(-0.56126422) q[1];
sx q[1];
rz(-0.47163481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4260226) q[0];
sx q[0];
rz(-1.7105118) q[0];
sx q[0];
rz(2.2989397) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.083375562) q[2];
sx q[2];
rz(-2.7365587) q[2];
sx q[2];
rz(2.0138559) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9593352) q[1];
sx q[1];
rz(-1.253009) q[1];
sx q[1];
rz(2.5885568) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0105693) q[3];
sx q[3];
rz(-2.4268055) q[3];
sx q[3];
rz(0.21372114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4857594) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(1.2539697) q[2];
rz(2.8433825) q[3];
sx q[3];
rz(-1.8658274) q[3];
sx q[3];
rz(-1.6037174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0373847) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(-2.3846159) q[0];
rz(1.0589927) q[1];
sx q[1];
rz(-1.6964361) q[1];
sx q[1];
rz(-0.34245488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9152074) q[0];
sx q[0];
rz(-1.3801551) q[0];
sx q[0];
rz(-1.9739499) q[0];
x q[1];
rz(-1.9497237) q[2];
sx q[2];
rz(-0.73816002) q[2];
sx q[2];
rz(-1.9110695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7360644) q[1];
sx q[1];
rz(-1.2031415) q[1];
sx q[1];
rz(2.6339732) q[1];
rz(-pi) q[2];
rz(0.60635546) q[3];
sx q[3];
rz(-0.49288756) q[3];
sx q[3];
rz(0.78745251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4464438) q[2];
sx q[2];
rz(-1.4651639) q[2];
sx q[2];
rz(-0.62015074) q[2];
rz(2.5946963) q[3];
sx q[3];
rz(-2.9345025) q[3];
sx q[3];
rz(-1.0774405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3243489) q[0];
sx q[0];
rz(-2.9381848) q[0];
sx q[0];
rz(1.2212344) q[0];
rz(-1.1034032) q[1];
sx q[1];
rz(-1.9788479) q[1];
sx q[1];
rz(2.3308636) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87830905) q[0];
sx q[0];
rz(-1.4118402) q[0];
sx q[0];
rz(-0.25860383) q[0];
x q[1];
rz(-2.2037494) q[2];
sx q[2];
rz(-2.1223017) q[2];
sx q[2];
rz(2.1126975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83119694) q[1];
sx q[1];
rz(-1.7094905) q[1];
sx q[1];
rz(3.0866361) q[1];
rz(-pi) q[2];
rz(2.8437988) q[3];
sx q[3];
rz(-1.9341365) q[3];
sx q[3];
rz(0.17223528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3409884) q[2];
sx q[2];
rz(-1.9438513) q[2];
sx q[2];
rz(1.9761337) q[2];
rz(-3.0321339) q[3];
sx q[3];
rz(-1.0215003) q[3];
sx q[3];
rz(-3.0808595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3133746) q[0];
sx q[0];
rz(-3.1309541) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(2.8984046) q[1];
sx q[1];
rz(-2.7311192) q[1];
sx q[1];
rz(2.5004255) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12578861) q[0];
sx q[0];
rz(-2.1282853) q[0];
sx q[0];
rz(2.1899738) q[0];
x q[1];
rz(-0.38491048) q[2];
sx q[2];
rz(-2.353047) q[2];
sx q[2];
rz(-2.6039518) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3873167) q[1];
sx q[1];
rz(-0.85016221) q[1];
sx q[1];
rz(-0.44447036) q[1];
rz(-pi) q[2];
rz(-0.23232122) q[3];
sx q[3];
rz(-1.182175) q[3];
sx q[3];
rz(0.53750932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0082561) q[2];
sx q[2];
rz(-1.9828321) q[2];
sx q[2];
rz(0.94375098) q[2];
rz(2.4933695) q[3];
sx q[3];
rz(-2.4002878) q[3];
sx q[3];
rz(-0.6664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42565313) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(-0.42914036) q[0];
rz(2.7375713) q[1];
sx q[1];
rz(-2.7262913) q[1];
sx q[1];
rz(-2.2228352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0149834) q[0];
sx q[0];
rz(-2.1169239) q[0];
sx q[0];
rz(1.2864497) q[0];
rz(2.0112923) q[2];
sx q[2];
rz(-0.28708946) q[2];
sx q[2];
rz(2.8237791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1240108) q[1];
sx q[1];
rz(-2.7582844) q[1];
sx q[1];
rz(0.95335754) q[1];
x q[2];
rz(-1.346896) q[3];
sx q[3];
rz(-0.32156518) q[3];
sx q[3];
rz(2.0617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.018844133) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(-1.3502236) q[2];
rz(-0.33251479) q[3];
sx q[3];
rz(-2.2612488) q[3];
sx q[3];
rz(-1.7727218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74828446) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(0.45676029) q[0];
rz(-1.4211753) q[1];
sx q[1];
rz(-1.8616118) q[1];
sx q[1];
rz(0.054904003) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0739912) q[0];
sx q[0];
rz(-1.5419863) q[0];
sx q[0];
rz(2.3037698) q[0];
rz(-1.1403528) q[2];
sx q[2];
rz(-1.9483742) q[2];
sx q[2];
rz(-2.0098639) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.6206606) q[1];
sx q[1];
rz(-2.0550613) q[1];
sx q[1];
rz(0.24543945) q[1];
rz(-pi) q[2];
rz(0.53590448) q[3];
sx q[3];
rz(-1.8632938) q[3];
sx q[3];
rz(-1.747594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2533337) q[2];
sx q[2];
rz(-2.4312225) q[2];
sx q[2];
rz(0.19676512) q[2];
rz(-1.0737859) q[3];
sx q[3];
rz(-1.1353227) q[3];
sx q[3];
rz(0.73608583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1573023) q[0];
sx q[0];
rz(-0.8534011) q[0];
sx q[0];
rz(-2.2421457) q[0];
rz(0.31117123) q[1];
sx q[1];
rz(-1.9322461) q[1];
sx q[1];
rz(1.7412294) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.78681) q[0];
sx q[0];
rz(-1.7799095) q[0];
sx q[0];
rz(2.1378921) q[0];
rz(-pi) q[1];
rz(2.0564931) q[2];
sx q[2];
rz(-1.3553095) q[2];
sx q[2];
rz(-2.4485916) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0728112) q[1];
sx q[1];
rz(-0.97807074) q[1];
sx q[1];
rz(2.0434421) q[1];
rz(2.4813249) q[3];
sx q[3];
rz(-0.83383027) q[3];
sx q[3];
rz(1.6136606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.060999) q[2];
sx q[2];
rz(-2.521535) q[2];
sx q[2];
rz(1.9798123) q[2];
rz(2.0628085) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(0.88206464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5494736) q[0];
sx q[0];
rz(-1.619512) q[0];
sx q[0];
rz(-1.6721538) q[0];
rz(0.1334162) q[1];
sx q[1];
rz(-1.7620371) q[1];
sx q[1];
rz(-0.98099991) q[1];
rz(-3.1405814) q[2];
sx q[2];
rz(-1.5173923) q[2];
sx q[2];
rz(-1.5569729) q[2];
rz(3.0117494) q[3];
sx q[3];
rz(-2.6363439) q[3];
sx q[3];
rz(1.5324203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

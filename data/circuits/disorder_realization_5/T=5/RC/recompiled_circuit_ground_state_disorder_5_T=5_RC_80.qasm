OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(-1.0125546) q[0];
sx q[0];
rz(-2.2192686) q[0];
rz(-0.84696472) q[1];
sx q[1];
rz(-1.6672517) q[1];
sx q[1];
rz(0.21811952) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6309752) q[0];
sx q[0];
rz(-1.6625064) q[0];
sx q[0];
rz(2.6789078) q[0];
rz(-pi) q[1];
rz(1.5297024) q[2];
sx q[2];
rz(-1.9485222) q[2];
sx q[2];
rz(-0.95313493) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0099758) q[1];
sx q[1];
rz(-2.9180315) q[1];
sx q[1];
rz(1.7625336) q[1];
rz(1.7401198) q[3];
sx q[3];
rz(-2.7062731) q[3];
sx q[3];
rz(3.0141413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0685136) q[2];
sx q[2];
rz(-1.928669) q[2];
sx q[2];
rz(-2.6050513) q[2];
rz(-1.0088629) q[3];
sx q[3];
rz(-2.3705685) q[3];
sx q[3];
rz(0.81396508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5040078) q[0];
sx q[0];
rz(-1.8300087) q[0];
sx q[0];
rz(-3.1245533) q[0];
rz(-1.0631961) q[1];
sx q[1];
rz(-2.3737962) q[1];
sx q[1];
rz(2.1868736) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9794036) q[0];
sx q[0];
rz(-1.5155751) q[0];
sx q[0];
rz(-2.9121141) q[0];
rz(1.8700908) q[2];
sx q[2];
rz(-1.7281282) q[2];
sx q[2];
rz(-0.89033876) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3054637) q[1];
sx q[1];
rz(-1.8082128) q[1];
sx q[1];
rz(1.3386471) q[1];
rz(-3.122284) q[3];
sx q[3];
rz(-1.5175765) q[3];
sx q[3];
rz(-1.9198708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9118328) q[2];
sx q[2];
rz(-0.91270295) q[2];
sx q[2];
rz(-2.0274053) q[2];
rz(1.1344502) q[3];
sx q[3];
rz(-0.19616923) q[3];
sx q[3];
rz(-0.1964868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5512307) q[0];
sx q[0];
rz(-0.95142618) q[0];
sx q[0];
rz(-0.1828585) q[0];
rz(0.081347801) q[1];
sx q[1];
rz(-0.75129879) q[1];
sx q[1];
rz(-0.61838165) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8066766) q[0];
sx q[0];
rz(-2.3551919) q[0];
sx q[0];
rz(1.4674835) q[0];
x q[1];
rz(-1.7612533) q[2];
sx q[2];
rz(-2.3229685) q[2];
sx q[2];
rz(0.10701767) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8733682) q[1];
sx q[1];
rz(-1.7098688) q[1];
sx q[1];
rz(-1.9130318) q[1];
x q[2];
rz(-0.86520536) q[3];
sx q[3];
rz(-0.83704797) q[3];
sx q[3];
rz(-3.1045543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0028093) q[2];
sx q[2];
rz(-1.323779) q[2];
sx q[2];
rz(2.7749824) q[2];
rz(1.2159411) q[3];
sx q[3];
rz(-1.9462908) q[3];
sx q[3];
rz(-2.478157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.65130305) q[0];
sx q[0];
rz(-1.2325352) q[0];
sx q[0];
rz(0.7269727) q[0];
rz(0.90826774) q[1];
sx q[1];
rz(-2.5636702) q[1];
sx q[1];
rz(1.04331) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9350727) q[0];
sx q[0];
rz(-1.5762202) q[0];
sx q[0];
rz(-0.66117735) q[0];
rz(-3.0798172) q[2];
sx q[2];
rz(-1.9456269) q[2];
sx q[2];
rz(1.0967364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2542782) q[1];
sx q[1];
rz(-1.320854) q[1];
sx q[1];
rz(2.0292086) q[1];
x q[2];
rz(2.7611012) q[3];
sx q[3];
rz(-2.0999319) q[3];
sx q[3];
rz(-0.38705119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72994453) q[2];
sx q[2];
rz(-2.8204155) q[2];
sx q[2];
rz(-1.8357065) q[2];
rz(0.11792396) q[3];
sx q[3];
rz(-0.4685466) q[3];
sx q[3];
rz(2.0809295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1082728) q[0];
sx q[0];
rz(-2.1554027) q[0];
sx q[0];
rz(1.6075851) q[0];
rz(-0.29574212) q[1];
sx q[1];
rz(-1.60138) q[1];
sx q[1];
rz(1.6827513) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.119232) q[0];
sx q[0];
rz(-1.5007449) q[0];
sx q[0];
rz(-0.82470597) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88250156) q[2];
sx q[2];
rz(-0.46445981) q[2];
sx q[2];
rz(-1.2683753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4319358) q[1];
sx q[1];
rz(-1.3742347) q[1];
sx q[1];
rz(-0.86039575) q[1];
x q[2];
rz(1.9146054) q[3];
sx q[3];
rz(-0.60887486) q[3];
sx q[3];
rz(0.46907779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7532588) q[2];
sx q[2];
rz(-0.39187852) q[2];
sx q[2];
rz(-0.64588109) q[2];
rz(-1.215747) q[3];
sx q[3];
rz(-1.682621) q[3];
sx q[3];
rz(-1.7997883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.129824) q[0];
sx q[0];
rz(-0.83368603) q[0];
sx q[0];
rz(-0.60761333) q[0];
rz(-1.1786849) q[1];
sx q[1];
rz(-1.8557529) q[1];
sx q[1];
rz(-2.7778621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4974198) q[0];
sx q[0];
rz(-1.8086642) q[0];
sx q[0];
rz(-2.1015443) q[0];
rz(1.2380139) q[2];
sx q[2];
rz(-1.6146891) q[2];
sx q[2];
rz(-1.3258758) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1235001) q[1];
sx q[1];
rz(-2.315633) q[1];
sx q[1];
rz(2.9781746) q[1];
rz(-pi) q[2];
rz(-0.97447864) q[3];
sx q[3];
rz(-1.5502366) q[3];
sx q[3];
rz(-2.1316776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58934775) q[2];
sx q[2];
rz(-3.037368) q[2];
sx q[2];
rz(-1.1927401) q[2];
rz(2.4097811) q[3];
sx q[3];
rz(-1.4251499) q[3];
sx q[3];
rz(1.4881136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8530387) q[0];
sx q[0];
rz(-2.9712501) q[0];
sx q[0];
rz(1.6188251) q[0];
rz(-1.4672) q[1];
sx q[1];
rz(-1.8312981) q[1];
sx q[1];
rz(-2.4640962) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554133) q[0];
sx q[0];
rz(-0.87550801) q[0];
sx q[0];
rz(-3.0762818) q[0];
rz(-pi) q[1];
rz(-0.67393731) q[2];
sx q[2];
rz(-0.29602414) q[2];
sx q[2];
rz(2.7715671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4851393) q[1];
sx q[1];
rz(-0.59018007) q[1];
sx q[1];
rz(0.82832576) q[1];
rz(1.7301637) q[3];
sx q[3];
rz(-2.2969118) q[3];
sx q[3];
rz(2.1350525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7447394) q[2];
sx q[2];
rz(-0.91529673) q[2];
sx q[2];
rz(-1.8358561) q[2];
rz(-0.33852494) q[3];
sx q[3];
rz(-1.1208231) q[3];
sx q[3];
rz(-2.4566076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6029538) q[0];
sx q[0];
rz(-2.9260577) q[0];
sx q[0];
rz(-0.18390528) q[0];
rz(-1.6869847) q[1];
sx q[1];
rz(-0.55211663) q[1];
sx q[1];
rz(-2.6457381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496807) q[0];
sx q[0];
rz(-0.56654585) q[0];
sx q[0];
rz(-2.0828155) q[0];
rz(-pi) q[1];
rz(0.74842986) q[2];
sx q[2];
rz(-2.1572621) q[2];
sx q[2];
rz(-0.85916729) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4019805) q[1];
sx q[1];
rz(-0.53369265) q[1];
sx q[1];
rz(2.3237565) q[1];
x q[2];
rz(3.0227376) q[3];
sx q[3];
rz(-2.1200075) q[3];
sx q[3];
rz(-2.1191747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0869861) q[2];
sx q[2];
rz(-0.84422529) q[2];
sx q[2];
rz(1.1620713) q[2];
rz(-1.6880796) q[3];
sx q[3];
rz(-0.96265692) q[3];
sx q[3];
rz(1.2801722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42089713) q[0];
sx q[0];
rz(-1.9891885) q[0];
sx q[0];
rz(2.5757117) q[0];
rz(2.7512918) q[1];
sx q[1];
rz(-1.5696328) q[1];
sx q[1];
rz(1.8338667) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5222675) q[0];
sx q[0];
rz(-0.59460708) q[0];
sx q[0];
rz(2.1451166) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0686223) q[2];
sx q[2];
rz(-0.42440571) q[2];
sx q[2];
rz(-0.54481943) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3257287) q[1];
sx q[1];
rz(-0.96885175) q[1];
sx q[1];
rz(1.25156) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2447137) q[3];
sx q[3];
rz(-0.19052902) q[3];
sx q[3];
rz(2.6485557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2719443) q[2];
sx q[2];
rz(-2.7248236) q[2];
sx q[2];
rz(2.1040253) q[2];
rz(-1.8218254) q[3];
sx q[3];
rz(-0.82873738) q[3];
sx q[3];
rz(-2.493609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8808402) q[0];
sx q[0];
rz(-0.97761959) q[0];
sx q[0];
rz(-2.4993437) q[0];
rz(1.3308659) q[1];
sx q[1];
rz(-0.86527491) q[1];
sx q[1];
rz(2.9790402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560682) q[0];
sx q[0];
rz(-2.3352288) q[0];
sx q[0];
rz(-2.1245405) q[0];
x q[1];
rz(-0.75679512) q[2];
sx q[2];
rz(-2.6399603) q[2];
sx q[2];
rz(-1.0197786) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7127258) q[1];
sx q[1];
rz(-1.3955294) q[1];
sx q[1];
rz(2.0574942) q[1];
rz(-pi) q[2];
rz(0.11998542) q[3];
sx q[3];
rz(-0.58661844) q[3];
sx q[3];
rz(-0.86931673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51005298) q[2];
sx q[2];
rz(-1.8659464) q[2];
sx q[2];
rz(-2.1007382) q[2];
rz(-2.1736274) q[3];
sx q[3];
rz(-2.0716045) q[3];
sx q[3];
rz(-2.4805243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3354202) q[0];
sx q[0];
rz(-1.5763043) q[0];
sx q[0];
rz(3.0446654) q[0];
rz(0.23282911) q[1];
sx q[1];
rz(-1.0284582) q[1];
sx q[1];
rz(-3.1340541) q[1];
rz(-2.3940968) q[2];
sx q[2];
rz(-1.1807673) q[2];
sx q[2];
rz(1.8058106) q[2];
rz(1.9540174) q[3];
sx q[3];
rz(-0.90771994) q[3];
sx q[3];
rz(-0.049605443) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

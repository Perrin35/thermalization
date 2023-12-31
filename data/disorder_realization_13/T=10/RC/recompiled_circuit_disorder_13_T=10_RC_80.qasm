OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(4.6306643) q[0];
sx q[0];
rz(10.319933) q[0];
rz(2.826638) q[1];
sx q[1];
rz(2.0575674) q[1];
sx q[1];
rz(10.881012) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2468949) q[0];
sx q[0];
rz(-0.1323192) q[0];
sx q[0];
rz(-1.7208862) q[0];
rz(-pi) q[1];
rz(2.7675682) q[2];
sx q[2];
rz(-2.0846539) q[2];
sx q[2];
rz(2.6467269) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40741062) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(2.461344) q[1];
x q[2];
rz(2.9772894) q[3];
sx q[3];
rz(-2.8176753) q[3];
sx q[3];
rz(1.2795554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3866117) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2125856) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(1.989495) q[0];
rz(1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5841056) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(1.4529865) q[0];
x q[1];
rz(2.9071964) q[2];
sx q[2];
rz(-2.1550551) q[2];
sx q[2];
rz(-0.61569475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86733782) q[1];
sx q[1];
rz(-1.6254289) q[1];
sx q[1];
rz(-1.0824624) q[1];
rz(-3.0456411) q[3];
sx q[3];
rz(-2.4437685) q[3];
sx q[3];
rz(-0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.058078893) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(0.2557959) q[2];
rz(1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26329041) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(0.27134744) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(0.43930611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8894316) q[0];
sx q[0];
rz(-1.541242) q[0];
sx q[0];
rz(1.6512524) q[0];
x q[1];
rz(-0.13621026) q[2];
sx q[2];
rz(-0.36362193) q[2];
sx q[2];
rz(1.9759535) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9648793) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(1.2769075) q[1];
rz(2.2112591) q[3];
sx q[3];
rz(-0.6593245) q[3];
sx q[3];
rz(1.7742771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1674041) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(2.9411194) q[2];
rz(-2.386507) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(-1.2444929) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(-2.2669852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4865206) q[0];
sx q[0];
rz(-0.97210303) q[0];
sx q[0];
rz(-0.18874164) q[0];
rz(0.34747296) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(-2.4328872) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7728459) q[1];
sx q[1];
rz(-2.6987942) q[1];
sx q[1];
rz(-2.5889791) q[1];
x q[2];
rz(1.4960257) q[3];
sx q[3];
rz(-1.6755591) q[3];
sx q[3];
rz(-1.9024224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-2.7757006) q[3];
sx q[3];
rz(1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(-2.5710035) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(-0.74329174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76737228) q[0];
sx q[0];
rz(-1.7905856) q[0];
sx q[0];
rz(0.90669294) q[0];
rz(-2.6402316) q[2];
sx q[2];
rz(-1.6622346) q[2];
sx q[2];
rz(-2.3102592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.78391) q[1];
sx q[1];
rz(-0.3158814) q[1];
sx q[1];
rz(2.3481579) q[1];
rz(-pi) q[2];
rz(2.8652142) q[3];
sx q[3];
rz(-2.31782) q[3];
sx q[3];
rz(0.36229047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9479998) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(2.8732079) q[2];
rz(-1.0466446) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(-2.823765) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.627581) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(2.6348689) q[0];
rz(0.2535893) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(1.0553029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21181606) q[0];
sx q[0];
rz(-1.7535216) q[0];
sx q[0];
rz(0.94434785) q[0];
x q[1];
rz(0.52168092) q[2];
sx q[2];
rz(-2.1573967) q[2];
sx q[2];
rz(-2.4751543) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.222059) q[1];
sx q[1];
rz(-2.2026081) q[1];
sx q[1];
rz(2.6421089) q[1];
rz(-pi) q[2];
rz(-0.25552337) q[3];
sx q[3];
rz(-2.3821085) q[3];
sx q[3];
rz(1.4042735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.48173299) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(2.4957538) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(-0.77254599) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(2.5678182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44874292) q[0];
sx q[0];
rz(-0.22073711) q[0];
sx q[0];
rz(-1.0462532) q[0];
rz(-pi) q[1];
rz(0.022106604) q[2];
sx q[2];
rz(-1.7709641) q[2];
sx q[2];
rz(-0.18891639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7635203) q[1];
sx q[1];
rz(-1.0135279) q[1];
sx q[1];
rz(0.65655638) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3366367) q[3];
sx q[3];
rz(-1.7171211) q[3];
sx q[3];
rz(2.9175266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8687826) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.1707206) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(1.8458813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1575748) q[0];
sx q[0];
rz(-2.5630953) q[0];
sx q[0];
rz(-0.63047854) q[0];
rz(-pi) q[1];
rz(-2.8700656) q[2];
sx q[2];
rz(-1.9121133) q[2];
sx q[2];
rz(1.5580387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.46591972) q[1];
sx q[1];
rz(-0.95341668) q[1];
sx q[1];
rz(-2.8663551) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0115764) q[3];
sx q[3];
rz(-1.5960346) q[3];
sx q[3];
rz(-2.4079635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7065113) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(0.56813017) q[2];
rz(0.98012296) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66529626) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(-2.4639159) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(2.303404) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8294551) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(1.8118993) q[0];
x q[1];
rz(-2.3269862) q[2];
sx q[2];
rz(-1.5216773) q[2];
sx q[2];
rz(-2.5533822) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9983478) q[1];
sx q[1];
rz(-1.76941) q[1];
sx q[1];
rz(2.3193588) q[1];
rz(1.4872876) q[3];
sx q[3];
rz(-1.828308) q[3];
sx q[3];
rz(-2.6058692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(1.8748803) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(-3.0468429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(-0.23751968) q[0];
rz(2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(0.231803) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7451413) q[0];
sx q[0];
rz(-1.6661577) q[0];
sx q[0];
rz(-2.0800637) q[0];
rz(0.78328697) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(0.17620262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69210359) q[1];
sx q[1];
rz(-2.8398501) q[1];
sx q[1];
rz(0.49441378) q[1];
rz(-pi) q[2];
rz(-1.8292571) q[3];
sx q[3];
rz(-0.70445326) q[3];
sx q[3];
rz(2.7750912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-2.0533662) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647472) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(2.5122535) q[2];
sx q[2];
rz(-2.7468801) q[2];
sx q[2];
rz(2.3135452) q[2];
rz(-0.18795342) q[3];
sx q[3];
rz(-1.5516075) q[3];
sx q[3];
rz(2.8967378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

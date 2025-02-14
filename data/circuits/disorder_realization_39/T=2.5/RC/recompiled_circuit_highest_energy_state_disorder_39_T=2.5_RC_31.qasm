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
rz(-2.0853618) q[0];
sx q[0];
rz(-0.27529278) q[0];
sx q[0];
rz(-3.0866301) q[0];
rz(-2.3995402) q[1];
sx q[1];
rz(-0.11062515) q[1];
sx q[1];
rz(-0.7551809) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36409828) q[0];
sx q[0];
rz(-0.75569442) q[0];
sx q[0];
rz(-0.022581935) q[0];
x q[1];
rz(-2.4651421) q[2];
sx q[2];
rz(-1.1764112) q[2];
sx q[2];
rz(2.6188322) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8480336) q[1];
sx q[1];
rz(-1.5747442) q[1];
sx q[1];
rz(-2.001954) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59261404) q[3];
sx q[3];
rz(-1.7120769) q[3];
sx q[3];
rz(-1.1326552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5521969) q[2];
sx q[2];
rz(-1.155921) q[2];
sx q[2];
rz(-1.7414306) q[2];
rz(-2.8081812) q[3];
sx q[3];
rz(-0.64781487) q[3];
sx q[3];
rz(2.6441372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19631504) q[0];
sx q[0];
rz(-0.59756398) q[0];
sx q[0];
rz(3.0481098) q[0];
rz(0.64158332) q[1];
sx q[1];
rz(-0.43951324) q[1];
sx q[1];
rz(-0.19215597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452602) q[0];
sx q[0];
rz(-2.30034) q[0];
sx q[0];
rz(1.6044046) q[0];
rz(2.9277705) q[2];
sx q[2];
rz(-1.7750083) q[2];
sx q[2];
rz(-1.363905) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9201492) q[1];
sx q[1];
rz(-2.0196947) q[1];
sx q[1];
rz(2.4549828) q[1];
x q[2];
rz(-0.3239969) q[3];
sx q[3];
rz(-1.5662875) q[3];
sx q[3];
rz(1.5229932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7265892) q[2];
sx q[2];
rz(-0.64565349) q[2];
sx q[2];
rz(2.1356706) q[2];
rz(3.1065324) q[3];
sx q[3];
rz(-0.56060767) q[3];
sx q[3];
rz(-2.3817271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1890892) q[0];
sx q[0];
rz(-1.3409505) q[0];
sx q[0];
rz(2.9553318) q[0];
rz(-2.8798036) q[1];
sx q[1];
rz(-1.2061385) q[1];
sx q[1];
rz(-2.632205) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68655339) q[0];
sx q[0];
rz(-1.5709366) q[0];
sx q[0];
rz(-1.554398) q[0];
rz(-pi) q[1];
rz(2.248412) q[2];
sx q[2];
rz(-1.2959157) q[2];
sx q[2];
rz(0.49472294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.81918908) q[1];
sx q[1];
rz(-1.6734646) q[1];
sx q[1];
rz(2.5346941) q[1];
rz(2.4548762) q[3];
sx q[3];
rz(-0.85609183) q[3];
sx q[3];
rz(1.2031632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.95106) q[2];
sx q[2];
rz(-2.4673927) q[2];
sx q[2];
rz(-0.29779693) q[2];
rz(-2.8262302) q[3];
sx q[3];
rz(-0.30839977) q[3];
sx q[3];
rz(1.4006008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0739338) q[0];
sx q[0];
rz(-0.44026259) q[0];
sx q[0];
rz(0.0010781188) q[0];
rz(0.06238097) q[1];
sx q[1];
rz(-1.9900813) q[1];
sx q[1];
rz(-1.6463702) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0001711) q[0];
sx q[0];
rz(-2.6323491) q[0];
sx q[0];
rz(-1.4285136) q[0];
rz(2.5220715) q[2];
sx q[2];
rz(-0.70370251) q[2];
sx q[2];
rz(1.0985398) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26186681) q[1];
sx q[1];
rz(-0.84964439) q[1];
sx q[1];
rz(2.7489064) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2132164) q[3];
sx q[3];
rz(-1.4331665) q[3];
sx q[3];
rz(-2.9127667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.262114) q[2];
sx q[2];
rz(-0.72747362) q[2];
sx q[2];
rz(1.0564085) q[2];
rz(-0.94111717) q[3];
sx q[3];
rz(-1.2818776) q[3];
sx q[3];
rz(3.0637872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5164186) q[0];
sx q[0];
rz(-2.43483) q[0];
sx q[0];
rz(-1.9551552) q[0];
rz(-2.2025853) q[1];
sx q[1];
rz(-0.72827315) q[1];
sx q[1];
rz(-0.77566385) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1224027) q[0];
sx q[0];
rz(-1.3870326) q[0];
sx q[0];
rz(1.8792436) q[0];
rz(-pi) q[1];
rz(-1.0177729) q[2];
sx q[2];
rz(-1.9423747) q[2];
sx q[2];
rz(-0.081173912) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.53117786) q[1];
sx q[1];
rz(-1.1503596) q[1];
sx q[1];
rz(-2.0095445) q[1];
rz(-0.19833447) q[3];
sx q[3];
rz(-2.9171067) q[3];
sx q[3];
rz(2.5189704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.080807216) q[2];
sx q[2];
rz(-2.5390517) q[2];
sx q[2];
rz(1.6281283) q[2];
rz(-1.3249409) q[3];
sx q[3];
rz(-0.63776582) q[3];
sx q[3];
rz(0.12346867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047664646) q[0];
sx q[0];
rz(-2.1340738) q[0];
sx q[0];
rz(-0.63543332) q[0];
rz(2.0280929) q[1];
sx q[1];
rz(-2.7058388) q[1];
sx q[1];
rz(-2.7472072) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52288234) q[0];
sx q[0];
rz(-0.11477192) q[0];
sx q[0];
rz(2.4023433) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7761781) q[2];
sx q[2];
rz(-0.49806777) q[2];
sx q[2];
rz(-0.77332449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.043076154) q[1];
sx q[1];
rz(-1.9253007) q[1];
sx q[1];
rz(0.47948821) q[1];
rz(-pi) q[2];
rz(-2.4784869) q[3];
sx q[3];
rz(-1.4713471) q[3];
sx q[3];
rz(0.15583043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8013132) q[2];
sx q[2];
rz(-1.3542513) q[2];
sx q[2];
rz(2.895943) q[2];
rz(0.73686016) q[3];
sx q[3];
rz(-2.1818678) q[3];
sx q[3];
rz(0.59011841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-3.0737632) q[0];
sx q[0];
rz(-1.2968061) q[0];
sx q[0];
rz(0.8514362) q[0];
rz(-0.40359452) q[1];
sx q[1];
rz(-2.5337296) q[1];
sx q[1];
rz(-3.0013989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4969585) q[0];
sx q[0];
rz(-1.2063364) q[0];
sx q[0];
rz(2.9011527) q[0];
x q[1];
rz(-1.6824016) q[2];
sx q[2];
rz(-1.3763469) q[2];
sx q[2];
rz(0.60332509) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.278001) q[1];
sx q[1];
rz(-0.74713445) q[1];
sx q[1];
rz(2.0160344) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1153489) q[3];
sx q[3];
rz(-1.035691) q[3];
sx q[3];
rz(-0.29087152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.95789528) q[2];
sx q[2];
rz(-0.87527466) q[2];
sx q[2];
rz(-2.9389006) q[2];
rz(1.1787339) q[3];
sx q[3];
rz(-0.26114994) q[3];
sx q[3];
rz(2.0451666) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75100768) q[0];
sx q[0];
rz(-3.106116) q[0];
sx q[0];
rz(1.9417199) q[0];
rz(-1.2212782) q[1];
sx q[1];
rz(-0.55322909) q[1];
sx q[1];
rz(1.663307) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8550028) q[0];
sx q[0];
rz(-0.012146771) q[0];
sx q[0];
rz(2.7783708) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7622819) q[2];
sx q[2];
rz(-2.0294445) q[2];
sx q[2];
rz(-0.5543602) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85046613) q[1];
sx q[1];
rz(-1.0134103) q[1];
sx q[1];
rz(-1.7859573) q[1];
x q[2];
rz(-2.6119328) q[3];
sx q[3];
rz(-0.33063355) q[3];
sx q[3];
rz(2.9042888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.025909802) q[2];
sx q[2];
rz(-2.2175711) q[2];
sx q[2];
rz(0.26205087) q[2];
rz(-0.15484364) q[3];
sx q[3];
rz(-2.9377929) q[3];
sx q[3];
rz(-3.1018597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994951) q[0];
sx q[0];
rz(-2.4423548) q[0];
sx q[0];
rz(-2.5788838) q[0];
rz(1.5744677) q[1];
sx q[1];
rz(-1.8561615) q[1];
sx q[1];
rz(-2.3318416) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8635559) q[0];
sx q[0];
rz(-0.42621729) q[0];
sx q[0];
rz(1.0025729) q[0];
x q[1];
rz(2.2024683) q[2];
sx q[2];
rz(-2.12471) q[2];
sx q[2];
rz(-2.6785452) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6470064) q[1];
sx q[1];
rz(-1.9368163) q[1];
sx q[1];
rz(-0.33509071) q[1];
rz(-3.0122981) q[3];
sx q[3];
rz(-2.2164549) q[3];
sx q[3];
rz(2.4753942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3385382) q[2];
sx q[2];
rz(-2.7092317) q[2];
sx q[2];
rz(1.7301249) q[2];
rz(2.9825397) q[3];
sx q[3];
rz(-1.8820857) q[3];
sx q[3];
rz(2.4771396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14801046) q[0];
sx q[0];
rz(-0.2033041) q[0];
sx q[0];
rz(2.5116442) q[0];
rz(2.3155164) q[1];
sx q[1];
rz(-0.45336205) q[1];
sx q[1];
rz(-2.1410543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8591908) q[0];
sx q[0];
rz(-1.0288718) q[0];
sx q[0];
rz(1.6247092) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8597322) q[2];
sx q[2];
rz(-2.0183767) q[2];
sx q[2];
rz(2.173169) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.40321585) q[1];
sx q[1];
rz(-2.1959702) q[1];
sx q[1];
rz(-1.0528864) q[1];
rz(-pi) q[2];
rz(-2.431683) q[3];
sx q[3];
rz(-1.1930575) q[3];
sx q[3];
rz(2.0722598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.91653812) q[2];
sx q[2];
rz(-0.30320898) q[2];
sx q[2];
rz(-2.0812601) q[2];
rz(0.66925085) q[3];
sx q[3];
rz(-2.4762912) q[3];
sx q[3];
rz(0.76702142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52436787) q[0];
sx q[0];
rz(-1.1001294) q[0];
sx q[0];
rz(-0.97468162) q[0];
rz(-2.2131447) q[1];
sx q[1];
rz(-1.7105553) q[1];
sx q[1];
rz(-1.1591844) q[1];
rz(-2.5090948) q[2];
sx q[2];
rz(-2.1618705) q[2];
sx q[2];
rz(2.4555997) q[2];
rz(-0.73406778) q[3];
sx q[3];
rz(-1.7597431) q[3];
sx q[3];
rz(-0.84398212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.2831777) q[0];
sx q[0];
rz(-1.5812961) q[0];
sx q[0];
rz(2.4726782) q[0];
rz(0.34612292) q[1];
sx q[1];
rz(-2.3200413) q[1];
sx q[1];
rz(-0.93924826) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6990841) q[0];
sx q[0];
rz(-0.48775771) q[0];
sx q[0];
rz(2.2217624) q[0];
x q[1];
rz(-1.6799567) q[2];
sx q[2];
rz(-2.7353035) q[2];
sx q[2];
rz(-0.10904212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2841442) q[1];
sx q[1];
rz(-1.836187) q[1];
sx q[1];
rz(-1.9678724) q[1];
rz(-0.40832728) q[3];
sx q[3];
rz(-0.75330594) q[3];
sx q[3];
rz(-1.7024794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5657438) q[2];
sx q[2];
rz(-2.1243024) q[2];
sx q[2];
rz(-3.1080833) q[2];
rz(1.2014028) q[3];
sx q[3];
rz(-1.7824495) q[3];
sx q[3];
rz(-2.0638154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1101542) q[0];
sx q[0];
rz(-0.25122508) q[0];
sx q[0];
rz(2.187619) q[0];
rz(-1.4454449) q[1];
sx q[1];
rz(-1.0547538) q[1];
sx q[1];
rz(1.9704069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8206089) q[0];
sx q[0];
rz(-2.1760328) q[0];
sx q[0];
rz(0.73802276) q[0];
rz(0.38925217) q[2];
sx q[2];
rz(-0.67485038) q[2];
sx q[2];
rz(-0.16216001) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9521683) q[1];
sx q[1];
rz(-1.599424) q[1];
sx q[1];
rz(2.9131123) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8071142) q[3];
sx q[3];
rz(-1.1568767) q[3];
sx q[3];
rz(-2.778307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90616068) q[2];
sx q[2];
rz(-1.7110598) q[2];
sx q[2];
rz(1.99235) q[2];
rz(-3.1084642) q[3];
sx q[3];
rz(-1.5002316) q[3];
sx q[3];
rz(2.5455425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0843435) q[0];
sx q[0];
rz(-0.79279041) q[0];
sx q[0];
rz(2.1113915) q[0];
rz(-1.239981) q[1];
sx q[1];
rz(-1.3571309) q[1];
sx q[1];
rz(2.1485567) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5035) q[0];
sx q[0];
rz(-2.6836392) q[0];
sx q[0];
rz(3.0328337) q[0];
x q[1];
rz(-2.9910261) q[2];
sx q[2];
rz(-0.67306256) q[2];
sx q[2];
rz(1.3526431) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5532304) q[1];
sx q[1];
rz(-0.85236406) q[1];
sx q[1];
rz(-2.3521982) q[1];
rz(-2.7168324) q[3];
sx q[3];
rz(-2.3094607) q[3];
sx q[3];
rz(0.55603851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9242636) q[2];
sx q[2];
rz(-2.0254878) q[2];
sx q[2];
rz(-0.35476157) q[2];
rz(-0.99700704) q[3];
sx q[3];
rz(-1.487178) q[3];
sx q[3];
rz(-0.94064373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(0.16911258) q[0];
sx q[0];
rz(-1.7153772) q[0];
sx q[0];
rz(1.1068363) q[0];
rz(-2.2726982) q[1];
sx q[1];
rz(-2.6241701) q[1];
sx q[1];
rz(-2.4443764) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1887564) q[0];
sx q[0];
rz(-1.7696912) q[0];
sx q[0];
rz(0.56629487) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4195693) q[2];
sx q[2];
rz(-1.486384) q[2];
sx q[2];
rz(-1.9868324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2576506) q[1];
sx q[1];
rz(-1.202794) q[1];
sx q[1];
rz(-1.9872527) q[1];
rz(-pi) q[2];
rz(-1.5460062) q[3];
sx q[3];
rz(-1.5953334) q[3];
sx q[3];
rz(2.8061574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8640459) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(-0.29402688) q[2];
rz(2.0781519) q[3];
sx q[3];
rz(-1.4592417) q[3];
sx q[3];
rz(-0.74240509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.71642891) q[0];
sx q[0];
rz(-2.9603781) q[0];
sx q[0];
rz(2.678405) q[0];
rz(2.7376392) q[1];
sx q[1];
rz(-1.5084167) q[1];
sx q[1];
rz(0.24872669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2866658) q[0];
sx q[0];
rz(-1.5724941) q[0];
sx q[0];
rz(1.3732852) q[0];
rz(-pi) q[1];
rz(-2.727365) q[2];
sx q[2];
rz(-2.5606206) q[2];
sx q[2];
rz(1.5609978) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.675093) q[1];
sx q[1];
rz(-2.6394301) q[1];
sx q[1];
rz(2.2999022) q[1];
x q[2];
rz(0.90051265) q[3];
sx q[3];
rz(-2.1888362) q[3];
sx q[3];
rz(-1.4423808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.8467466) q[2];
sx q[2];
rz(-1.8241901) q[2];
sx q[2];
rz(1.3999636) q[2];
rz(1.8585662) q[3];
sx q[3];
rz(-1.7581519) q[3];
sx q[3];
rz(2.9156901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95752174) q[0];
sx q[0];
rz(-2.3029843) q[0];
sx q[0];
rz(-0.34570178) q[0];
rz(-0.050994571) q[1];
sx q[1];
rz(-1.9146999) q[1];
sx q[1];
rz(-0.7720224) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22020082) q[0];
sx q[0];
rz(-1.8225551) q[0];
sx q[0];
rz(2.779106) q[0];
rz(-2.2228681) q[2];
sx q[2];
rz(-1.6479392) q[2];
sx q[2];
rz(-0.89706286) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.009338) q[1];
sx q[1];
rz(-1.8937832) q[1];
sx q[1];
rz(-0.96829523) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0984678) q[3];
sx q[3];
rz(-1.3195795) q[3];
sx q[3];
rz(-0.71300292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2766075) q[2];
sx q[2];
rz(-1.425068) q[2];
sx q[2];
rz(2.9841606) q[2];
rz(-2.7031247) q[3];
sx q[3];
rz(-0.74123588) q[3];
sx q[3];
rz(0.95611519) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081414374) q[0];
sx q[0];
rz(-2.0745451) q[0];
sx q[0];
rz(9/(11*pi)) q[0];
rz(2.5632437) q[1];
sx q[1];
rz(-1.7702421) q[1];
sx q[1];
rz(2.9885805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4082913) q[0];
sx q[0];
rz(-0.75364764) q[0];
sx q[0];
rz(-0.80259364) q[0];
x q[1];
rz(-3.0566759) q[2];
sx q[2];
rz(-1.474829) q[2];
sx q[2];
rz(1.7964448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2305829) q[1];
sx q[1];
rz(-2.1248159) q[1];
sx q[1];
rz(-1.0826712) q[1];
rz(-pi) q[2];
rz(2.8197391) q[3];
sx q[3];
rz(-2.0087089) q[3];
sx q[3];
rz(1.9066325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.962062) q[2];
sx q[2];
rz(-0.099055812) q[2];
sx q[2];
rz(2.8774101) q[2];
rz(-1.8386486) q[3];
sx q[3];
rz(-1.3774201) q[3];
sx q[3];
rz(-2.0343659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95689479) q[0];
sx q[0];
rz(-2.1883924) q[0];
sx q[0];
rz(-1.4554998) q[0];
rz(-0.9616583) q[1];
sx q[1];
rz(-1.7627629) q[1];
sx q[1];
rz(0.89967322) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6379999) q[0];
sx q[0];
rz(-1.3595306) q[0];
sx q[0];
rz(0.46163033) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3191965) q[2];
sx q[2];
rz(-2.4008022) q[2];
sx q[2];
rz(-1.0902001) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74994319) q[1];
sx q[1];
rz(-0.98403059) q[1];
sx q[1];
rz(-3.090078) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7387271) q[3];
sx q[3];
rz(-0.21997082) q[3];
sx q[3];
rz(-2.0988362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4453033) q[2];
sx q[2];
rz(-0.60089198) q[2];
sx q[2];
rz(0.75330934) q[2];
rz(2.8478012) q[3];
sx q[3];
rz(-1.1283504) q[3];
sx q[3];
rz(-2.0297348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488572) q[0];
sx q[0];
rz(-0.98965544) q[0];
sx q[0];
rz(-2.2897172) q[0];
rz(1.4082255) q[1];
sx q[1];
rz(-1.4829166) q[1];
sx q[1];
rz(2.7588989) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9155884) q[0];
sx q[0];
rz(-1.3289641) q[0];
sx q[0];
rz(0.42892021) q[0];
x q[1];
rz(2.9054602) q[2];
sx q[2];
rz(-2.3575767) q[2];
sx q[2];
rz(1.2003984) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1074658) q[1];
sx q[1];
rz(-1.6686979) q[1];
sx q[1];
rz(0.30219725) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9858573) q[3];
sx q[3];
rz(-0.49113516) q[3];
sx q[3];
rz(-2.8409426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3771628) q[2];
sx q[2];
rz(-0.96401507) q[2];
sx q[2];
rz(0.60297472) q[2];
rz(2.073334) q[3];
sx q[3];
rz(-1.7030741) q[3];
sx q[3];
rz(-2.300613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0076440796) q[0];
sx q[0];
rz(-1.4643865) q[0];
sx q[0];
rz(-0.35368791) q[0];
rz(1.1324646) q[1];
sx q[1];
rz(-0.9681038) q[1];
sx q[1];
rz(-2.7521334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2089473) q[0];
sx q[0];
rz(-1.07558) q[0];
sx q[0];
rz(1.8344648) q[0];
rz(-pi) q[1];
rz(-3.0385397) q[2];
sx q[2];
rz(-1.3259058) q[2];
sx q[2];
rz(2.6239397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8824401) q[1];
sx q[1];
rz(-2.5969567) q[1];
sx q[1];
rz(-0.55014054) q[1];
rz(-pi) q[2];
rz(1.6234446) q[3];
sx q[3];
rz(-2.0984167) q[3];
sx q[3];
rz(1.2155346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0103717) q[2];
sx q[2];
rz(-1.5067357) q[2];
sx q[2];
rz(-1.9949251) q[2];
rz(1.5366588) q[3];
sx q[3];
rz(-0.48303548) q[3];
sx q[3];
rz(-2.7893132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1491886) q[0];
sx q[0];
rz(-0.70269361) q[0];
sx q[0];
rz(0.50444034) q[0];
rz(-3.0814677) q[1];
sx q[1];
rz(-2.6838214) q[1];
sx q[1];
rz(2.6837742) q[1];
rz(2.1679466) q[2];
sx q[2];
rz(-0.13046593) q[2];
sx q[2];
rz(2.2400357) q[2];
rz(-2.2452284) q[3];
sx q[3];
rz(-1.5980362) q[3];
sx q[3];
rz(-2.7841795) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

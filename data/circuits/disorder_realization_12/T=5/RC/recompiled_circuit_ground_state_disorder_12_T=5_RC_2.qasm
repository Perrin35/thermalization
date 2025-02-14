OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0491068) q[0];
sx q[0];
rz(-1.8422814) q[0];
sx q[0];
rz(0.045825034) q[0];
rz(1.2996281) q[1];
sx q[1];
rz(-1.8027432) q[1];
sx q[1];
rz(1.2531228) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9535256) q[0];
sx q[0];
rz(-1.3261516) q[0];
sx q[0];
rz(1.1172386) q[0];
rz(-2.9304977) q[2];
sx q[2];
rz(-0.20583868) q[2];
sx q[2];
rz(-1.7072276) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8032082) q[1];
sx q[1];
rz(-2.3178891) q[1];
sx q[1];
rz(-1.4080774) q[1];
rz(-pi) q[2];
rz(1.5398938) q[3];
sx q[3];
rz(-3.0140844) q[3];
sx q[3];
rz(2.7650583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32156285) q[2];
sx q[2];
rz(-1.4048615) q[2];
sx q[2];
rz(-0.31537867) q[2];
rz(3.0553715) q[3];
sx q[3];
rz(-0.092970522) q[3];
sx q[3];
rz(0.84403795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.740199) q[0];
sx q[0];
rz(-1.4305776) q[0];
sx q[0];
rz(2.8061818) q[0];
rz(-2.7566578) q[1];
sx q[1];
rz(-2.3454869) q[1];
sx q[1];
rz(1.2892494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3649454) q[0];
sx q[0];
rz(-1.9595265) q[0];
sx q[0];
rz(-1.5035648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.094893242) q[2];
sx q[2];
rz(-0.93145639) q[2];
sx q[2];
rz(-2.4415093) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8131539) q[1];
sx q[1];
rz(-0.48373172) q[1];
sx q[1];
rz(0.64240047) q[1];
x q[2];
rz(0.13522526) q[3];
sx q[3];
rz(-1.3298508) q[3];
sx q[3];
rz(1.2775482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59911597) q[2];
sx q[2];
rz(-1.7581538) q[2];
sx q[2];
rz(-1.6790338) q[2];
rz(2.2595432) q[3];
sx q[3];
rz(-2.4817011) q[3];
sx q[3];
rz(1.6409469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7733234) q[0];
sx q[0];
rz(-2.6550846) q[0];
sx q[0];
rz(-0.64273709) q[0];
rz(2.2679988) q[1];
sx q[1];
rz(-1.3106376) q[1];
sx q[1];
rz(0.98683039) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1571647) q[0];
sx q[0];
rz(-1.4207977) q[0];
sx q[0];
rz(0.87777975) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6880694) q[2];
sx q[2];
rz(-0.48732317) q[2];
sx q[2];
rz(-0.27299689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22443188) q[1];
sx q[1];
rz(-2.3244435) q[1];
sx q[1];
rz(-2.3520873) q[1];
rz(-pi) q[2];
rz(1.1961206) q[3];
sx q[3];
rz(-1.6171394) q[3];
sx q[3];
rz(-2.1425193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4387536) q[2];
sx q[2];
rz(-0.025391014) q[2];
sx q[2];
rz(2.0849483) q[2];
rz(-1.7993401) q[3];
sx q[3];
rz(-1.4547576) q[3];
sx q[3];
rz(-2.2630283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037647) q[0];
sx q[0];
rz(-0.28452888) q[0];
sx q[0];
rz(1.9985265) q[0];
rz(2.4567538) q[1];
sx q[1];
rz(-1.3416483) q[1];
sx q[1];
rz(0.24982223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78501084) q[0];
sx q[0];
rz(-2.2306721) q[0];
sx q[0];
rz(-0.59156811) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1353829) q[2];
sx q[2];
rz(-2.5673625) q[2];
sx q[2];
rz(-1.8956313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.494809) q[1];
sx q[1];
rz(-0.80776359) q[1];
sx q[1];
rz(-2.450472) q[1];
rz(-3.1400348) q[3];
sx q[3];
rz(-1.8127725) q[3];
sx q[3];
rz(-0.93780692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7217241) q[2];
sx q[2];
rz(-1.7896174) q[2];
sx q[2];
rz(1.8458337) q[2];
rz(-1.3755679) q[3];
sx q[3];
rz(-1.7121168) q[3];
sx q[3];
rz(-0.099055722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-2.7527723) q[0];
sx q[0];
rz(-1.3714014) q[0];
sx q[0];
rz(-0.8465299) q[0];
rz(-1.9580152) q[1];
sx q[1];
rz(-1.9639683) q[1];
sx q[1];
rz(1.902098) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0087933) q[0];
sx q[0];
rz(-1.4504878) q[0];
sx q[0];
rz(1.8254721) q[0];
x q[1];
rz(-0.47649224) q[2];
sx q[2];
rz(-2.7600754) q[2];
sx q[2];
rz(-3.1135349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.848865) q[1];
sx q[1];
rz(-1.5588817) q[1];
sx q[1];
rz(1.6267651) q[1];
rz(1.6587448) q[3];
sx q[3];
rz(-0.59895688) q[3];
sx q[3];
rz(-0.14725895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6401297) q[2];
sx q[2];
rz(-1.5093466) q[2];
sx q[2];
rz(-0.35422361) q[2];
rz(0.7192449) q[3];
sx q[3];
rz(-0.57643276) q[3];
sx q[3];
rz(-0.80938068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3767913) q[0];
sx q[0];
rz(-0.23623315) q[0];
sx q[0];
rz(0.34673044) q[0];
rz(-0.72225839) q[1];
sx q[1];
rz(-1.6112695) q[1];
sx q[1];
rz(-2.9647656) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5895329) q[0];
sx q[0];
rz(-1.1101652) q[0];
sx q[0];
rz(0.15710196) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7406101) q[2];
sx q[2];
rz(-0.50504518) q[2];
sx q[2];
rz(1.9521879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4011098) q[1];
sx q[1];
rz(-1.4663457) q[1];
sx q[1];
rz(-2.9354276) q[1];
rz(-0.10156722) q[3];
sx q[3];
rz(-2.1312765) q[3];
sx q[3];
rz(1.5965726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5727545) q[2];
sx q[2];
rz(-0.73358959) q[2];
sx q[2];
rz(-1.0432165) q[2];
rz(-0.15626945) q[3];
sx q[3];
rz(-1.0633435) q[3];
sx q[3];
rz(0.094303057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909376) q[0];
sx q[0];
rz(-2.1112878) q[0];
sx q[0];
rz(-1.913273) q[0];
rz(-0.43371513) q[1];
sx q[1];
rz(-2.3408196) q[1];
sx q[1];
rz(2.0965651) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50253326) q[0];
sx q[0];
rz(-2.9307588) q[0];
sx q[0];
rz(1.2913713) q[0];
x q[1];
rz(0.47335923) q[2];
sx q[2];
rz(-0.88027647) q[2];
sx q[2];
rz(2.8831633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6562798) q[1];
sx q[1];
rz(-1.4569786) q[1];
sx q[1];
rz(-0.95358221) q[1];
rz(-0.52900903) q[3];
sx q[3];
rz(-0.85011357) q[3];
sx q[3];
rz(-2.1095654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23994437) q[2];
sx q[2];
rz(-2.1176691) q[2];
sx q[2];
rz(-0.18590064) q[2];
rz(-1.2402041) q[3];
sx q[3];
rz(-1.600772) q[3];
sx q[3];
rz(1.0143636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1321201) q[0];
sx q[0];
rz(-1.2484231) q[0];
sx q[0];
rz(2.9820719) q[0];
rz(0.80687833) q[1];
sx q[1];
rz(-1.9069549) q[1];
sx q[1];
rz(0.0085011403) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59183622) q[0];
sx q[0];
rz(-1.7248132) q[0];
sx q[0];
rz(1.9817673) q[0];
rz(1.464017) q[2];
sx q[2];
rz(-0.77645436) q[2];
sx q[2];
rz(-1.8242893) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2841683) q[1];
sx q[1];
rz(-1.4393974) q[1];
sx q[1];
rz(-2.9144276) q[1];
rz(2.9656677) q[3];
sx q[3];
rz(-2.0485544) q[3];
sx q[3];
rz(-2.0328409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8622417) q[2];
sx q[2];
rz(-1.3631577) q[2];
sx q[2];
rz(-0.19374338) q[2];
rz(1.5396317) q[3];
sx q[3];
rz(-1.964147) q[3];
sx q[3];
rz(-1.1295454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2368471) q[0];
sx q[0];
rz(-1.3847677) q[0];
sx q[0];
rz(-2.4054476) q[0];
rz(0.88788095) q[1];
sx q[1];
rz(-2.6887951) q[1];
sx q[1];
rz(-3.091541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5745478) q[0];
sx q[0];
rz(-2.7549681) q[0];
sx q[0];
rz(2.6599314) q[0];
rz(-2.8456837) q[2];
sx q[2];
rz(-2.8136133) q[2];
sx q[2];
rz(2.8519972) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5347595) q[1];
sx q[1];
rz(-1.2616565) q[1];
sx q[1];
rz(-1.3571795) q[1];
rz(2.1878236) q[3];
sx q[3];
rz(-2.0514384) q[3];
sx q[3];
rz(-2.1565089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2729518) q[2];
sx q[2];
rz(-1.7014039) q[2];
sx q[2];
rz(0.0090553332) q[2];
rz(2.791259) q[3];
sx q[3];
rz(-0.30600268) q[3];
sx q[3];
rz(2.8370324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7250799) q[0];
sx q[0];
rz(-1.059499) q[0];
sx q[0];
rz(-0.9730202) q[0];
rz(1.579772) q[1];
sx q[1];
rz(-0.90619722) q[1];
sx q[1];
rz(2.3320893) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79671395) q[0];
sx q[0];
rz(-2.2911706) q[0];
sx q[0];
rz(1.0864837) q[0];
rz(1.3415178) q[2];
sx q[2];
rz(-1.3827168) q[2];
sx q[2];
rz(0.089872472) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36724801) q[1];
sx q[1];
rz(-1.4092236) q[1];
sx q[1];
rz(-1.104924) q[1];
rz(0.46790917) q[3];
sx q[3];
rz(-1.0413326) q[3];
sx q[3];
rz(2.5264945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96106625) q[2];
sx q[2];
rz(-2.3193391) q[2];
sx q[2];
rz(0.52214885) q[2];
rz(2.7850049) q[3];
sx q[3];
rz(-0.57727376) q[3];
sx q[3];
rz(3.0780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4840354) q[0];
sx q[0];
rz(-0.35150305) q[0];
sx q[0];
rz(1.2765314) q[0];
rz(2.7784078) q[1];
sx q[1];
rz(-2.6138432) q[1];
sx q[1];
rz(1.6539727) q[1];
rz(2.364306) q[2];
sx q[2];
rz(-1.1697792) q[2];
sx q[2];
rz(-0.32464632) q[2];
rz(2.667712) q[3];
sx q[3];
rz(-0.6345748) q[3];
sx q[3];
rz(2.0848772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

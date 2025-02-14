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
rz(1.8013826) q[0];
sx q[0];
rz(-0.27685452) q[0];
sx q[0];
rz(2.1028331) q[0];
rz(-0.34173319) q[1];
sx q[1];
rz(-2.3050397) q[1];
sx q[1];
rz(0.41681448) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78535372) q[0];
sx q[0];
rz(-2.2921341) q[0];
sx q[0];
rz(-1.712404) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7567467) q[2];
sx q[2];
rz(-0.59078465) q[2];
sx q[2];
rz(-2.7123775) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8326679) q[1];
sx q[1];
rz(-1.0005669) q[1];
sx q[1];
rz(-0.69242386) q[1];
x q[2];
rz(-1.5384212) q[3];
sx q[3];
rz(-1.4093295) q[3];
sx q[3];
rz(2.3700455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9716399) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(-1.2539585) q[2];
rz(-1.5597255) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1935683) q[0];
sx q[0];
rz(-1.7727611) q[0];
sx q[0];
rz(0.76876202) q[0];
rz(-1.3668775) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(-1.0381402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4878142) q[0];
sx q[0];
rz(-0.013049203) q[0];
sx q[0];
rz(2.292146) q[0];
rz(2.4203796) q[2];
sx q[2];
rz(-1.8497457) q[2];
sx q[2];
rz(-2.9245289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8703192) q[1];
sx q[1];
rz(-2.4044577) q[1];
sx q[1];
rz(2.5802274) q[1];
rz(0.019006218) q[3];
sx q[3];
rz(-1.0822625) q[3];
sx q[3];
rz(1.8754539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64608964) q[2];
sx q[2];
rz(-1.2868519) q[2];
sx q[2];
rz(2.3213279) q[2];
rz(-0.25343728) q[3];
sx q[3];
rz(-0.35496747) q[3];
sx q[3];
rz(0.36111116) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8609817) q[0];
sx q[0];
rz(-2.6237223) q[0];
sx q[0];
rz(-2.2542727) q[0];
rz(0.53030983) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(-1.1393772) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054087669) q[0];
sx q[0];
rz(-1.2257135) q[0];
sx q[0];
rz(-1.7214543) q[0];
x q[1];
rz(-0.33311756) q[2];
sx q[2];
rz(-1.0651759) q[2];
sx q[2];
rz(-0.62084711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92147747) q[1];
sx q[1];
rz(-2.8954828) q[1];
sx q[1];
rz(3.0819478) q[1];
x q[2];
rz(0.51872571) q[3];
sx q[3];
rz(-1.5253807) q[3];
sx q[3];
rz(-2.1741427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.016971074) q[2];
sx q[2];
rz(-1.6861702) q[2];
sx q[2];
rz(-2.7023081) q[2];
rz(0.47075054) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(1.9973756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73862326) q[0];
sx q[0];
rz(-2.2082177) q[0];
sx q[0];
rz(3.0261107) q[0];
rz(0.11784095) q[1];
sx q[1];
rz(-1.9385447) q[1];
sx q[1];
rz(-1.3062564) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8354315) q[0];
sx q[0];
rz(-1.8580313) q[0];
sx q[0];
rz(-1.1348508) q[0];
rz(-0.47331402) q[2];
sx q[2];
rz(-0.71719786) q[2];
sx q[2];
rz(-2.0462478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.072767898) q[1];
sx q[1];
rz(-0.95736527) q[1];
sx q[1];
rz(-2.3062506) q[1];
rz(0.0751817) q[3];
sx q[3];
rz(-2.1422221) q[3];
sx q[3];
rz(0.21440766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8370886) q[2];
sx q[2];
rz(-1.3136761) q[2];
sx q[2];
rz(2.9869249) q[2];
rz(-1.6020487) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(-0.60607564) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57347572) q[0];
sx q[0];
rz(-1.0433759) q[0];
sx q[0];
rz(-0.83531761) q[0];
rz(-0.71594816) q[1];
sx q[1];
rz(-2.7138111) q[1];
sx q[1];
rz(-3.0221525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1560005) q[0];
sx q[0];
rz(-1.5876666) q[0];
sx q[0];
rz(1.5149679) q[0];
x q[1];
rz(-1.4177464) q[2];
sx q[2];
rz(-2.0156392) q[2];
sx q[2];
rz(-0.21013021) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2970256) q[1];
sx q[1];
rz(-0.78999619) q[1];
sx q[1];
rz(-2.4639735) q[1];
x q[2];
rz(-0.53262696) q[3];
sx q[3];
rz(-1.2029193) q[3];
sx q[3];
rz(1.0581349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9014088) q[2];
sx q[2];
rz(-0.26075026) q[2];
sx q[2];
rz(3.0805947) q[2];
rz(3.0933464) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(-0.72859305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.401684) q[0];
sx q[0];
rz(-2.4329199) q[0];
sx q[0];
rz(-1.1309062) q[0];
rz(-0.4785969) q[1];
sx q[1];
rz(-2.5391948) q[1];
sx q[1];
rz(1.7344249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3307307) q[0];
sx q[0];
rz(-1.5172345) q[0];
sx q[0];
rz(1.5723482) q[0];
rz(-0.84848225) q[2];
sx q[2];
rz(-1.7847848) q[2];
sx q[2];
rz(0.011836476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1651014) q[1];
sx q[1];
rz(-2.8860984) q[1];
sx q[1];
rz(2.0352023) q[1];
x q[2];
rz(2.4396517) q[3];
sx q[3];
rz(-0.27595584) q[3];
sx q[3];
rz(-2.4380033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5799134) q[2];
sx q[2];
rz(-0.40739569) q[2];
sx q[2];
rz(-1.178406) q[2];
rz(-2.0922349) q[3];
sx q[3];
rz(-0.5368084) q[3];
sx q[3];
rz(1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2105763) q[0];
sx q[0];
rz(-1.9518305) q[0];
sx q[0];
rz(1.3354906) q[0];
rz(1.8203075) q[1];
sx q[1];
rz(-2.0704465) q[1];
sx q[1];
rz(0.30805045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9481407) q[0];
sx q[0];
rz(-1.7823671) q[0];
sx q[0];
rz(-0.39724644) q[0];
rz(-pi) q[1];
rz(-0.51527889) q[2];
sx q[2];
rz(-1.9204428) q[2];
sx q[2];
rz(-0.51008979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96770331) q[1];
sx q[1];
rz(-2.1427665) q[1];
sx q[1];
rz(-2.553945) q[1];
rz(-1.1847897) q[3];
sx q[3];
rz(-1.526745) q[3];
sx q[3];
rz(3.0795627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41219741) q[2];
sx q[2];
rz(-0.9407548) q[2];
sx q[2];
rz(-3.0103053) q[2];
rz(-2.4042551) q[3];
sx q[3];
rz(-1.3056583) q[3];
sx q[3];
rz(-2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8312296) q[0];
sx q[0];
rz(-2.0518301) q[0];
sx q[0];
rz(-2.6112153) q[0];
rz(-1.4250379) q[1];
sx q[1];
rz(-2.0030463) q[1];
sx q[1];
rz(-2.4200965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8222161) q[0];
sx q[0];
rz(-2.0160897) q[0];
sx q[0];
rz(0.67586835) q[0];
rz(-pi) q[1];
rz(1.4230096) q[2];
sx q[2];
rz(-2.035454) q[2];
sx q[2];
rz(-2.8140659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87807314) q[1];
sx q[1];
rz(-1.7612447) q[1];
sx q[1];
rz(1.3506804) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3590478) q[3];
sx q[3];
rz(-1.492332) q[3];
sx q[3];
rz(-0.162938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43712744) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(-1.4777769) q[2];
rz(1.1635228) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32996938) q[0];
sx q[0];
rz(-1.7003308) q[0];
sx q[0];
rz(-0.60923088) q[0];
rz(1.5513264) q[1];
sx q[1];
rz(-2.8158999) q[1];
sx q[1];
rz(-1.4564266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.698899) q[0];
sx q[0];
rz(-2.6999268) q[0];
sx q[0];
rz(-0.95516686) q[0];
rz(-3.0266043) q[2];
sx q[2];
rz(-0.48795944) q[2];
sx q[2];
rz(1.8166627) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2472154) q[1];
sx q[1];
rz(-0.4900107) q[1];
sx q[1];
rz(-1.2213329) q[1];
rz(-pi) q[2];
rz(0.88202694) q[3];
sx q[3];
rz(-0.85245624) q[3];
sx q[3];
rz(0.73757987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41813254) q[2];
sx q[2];
rz(-2.2432566) q[2];
sx q[2];
rz(2.2612803) q[2];
rz(2.9355925) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(-2.6515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3708165) q[0];
sx q[0];
rz(-2.4801065) q[0];
sx q[0];
rz(0.01195512) q[0];
rz(-1.5097584) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(-2.5451122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76145455) q[0];
sx q[0];
rz(-1.5517457) q[0];
sx q[0];
rz(2.6721902) q[0];
x q[1];
rz(-0.61052236) q[2];
sx q[2];
rz(-2.7075665) q[2];
sx q[2];
rz(-2.4422925) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47190753) q[1];
sx q[1];
rz(-2.0457256) q[1];
sx q[1];
rz(0.088598786) q[1];
rz(-pi) q[2];
rz(2.8438472) q[3];
sx q[3];
rz(-1.7216873) q[3];
sx q[3];
rz(2.6835364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.163588) q[2];
sx q[2];
rz(-2.1808193) q[2];
sx q[2];
rz(2.9445924) q[2];
rz(1.152285) q[3];
sx q[3];
rz(-2.9983493) q[3];
sx q[3];
rz(0.44767374) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86354179) q[0];
sx q[0];
rz(-2.1320237) q[0];
sx q[0];
rz(-0.19436819) q[0];
rz(-0.20881431) q[1];
sx q[1];
rz(-1.5369692) q[1];
sx q[1];
rz(2.1280638) q[1];
rz(2.2603358) q[2];
sx q[2];
rz(-1.4980346) q[2];
sx q[2];
rz(-3.0537506) q[2];
rz(-2.6326451) q[3];
sx q[3];
rz(-1.043368) q[3];
sx q[3];
rz(2.1989078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

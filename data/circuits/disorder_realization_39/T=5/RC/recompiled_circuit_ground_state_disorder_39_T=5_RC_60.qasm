OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47867632) q[0];
sx q[0];
rz(7.8362099) q[0];
sx q[0];
rz(10.683164) q[0];
rz(1.2619184) q[1];
sx q[1];
rz(-2.6231397) q[1];
sx q[1];
rz(-0.060103091) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4661048) q[0];
sx q[0];
rz(-1.6610613) q[0];
sx q[0];
rz(1.6016927) q[0];
x q[1];
rz(2.8608039) q[2];
sx q[2];
rz(-0.5092237) q[2];
sx q[2];
rz(-1.4851242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93906885) q[1];
sx q[1];
rz(-1.9323404) q[1];
sx q[1];
rz(0.099379813) q[1];
rz(0.97385554) q[3];
sx q[3];
rz(-1.1789448) q[3];
sx q[3];
rz(2.4490956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95187676) q[2];
sx q[2];
rz(-2.3592301) q[2];
sx q[2];
rz(-2.961109) q[2];
rz(-2.8372724) q[3];
sx q[3];
rz(-2.2173827) q[3];
sx q[3];
rz(-0.2271823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75905281) q[0];
sx q[0];
rz(-0.57342425) q[0];
sx q[0];
rz(0.10391129) q[0];
rz(2.3567764) q[1];
sx q[1];
rz(-1.3094614) q[1];
sx q[1];
rz(-0.58194247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79051414) q[0];
sx q[0];
rz(-2.1026372) q[0];
sx q[0];
rz(1.5966589) q[0];
rz(-2.632276) q[2];
sx q[2];
rz(-1.2602206) q[2];
sx q[2];
rz(-0.79307014) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71680333) q[1];
sx q[1];
rz(-1.5967073) q[1];
sx q[1];
rz(2.0477363) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1725015) q[3];
sx q[3];
rz(-1.1886667) q[3];
sx q[3];
rz(-0.39234658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.71501422) q[2];
sx q[2];
rz(-0.38807401) q[2];
sx q[2];
rz(-1.0283872) q[2];
rz(0.55108023) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(-0.50857956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55117637) q[0];
sx q[0];
rz(-1.506378) q[0];
sx q[0];
rz(-1.2580385) q[0];
rz(-0.60351562) q[1];
sx q[1];
rz(-1.3110833) q[1];
sx q[1];
rz(-0.18361941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063832961) q[0];
sx q[0];
rz(-0.49657425) q[0];
sx q[0];
rz(-0.99061503) q[0];
x q[1];
rz(0.48086353) q[2];
sx q[2];
rz(-1.6893975) q[2];
sx q[2];
rz(3.1080217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.078122) q[1];
sx q[1];
rz(-0.87492157) q[1];
sx q[1];
rz(-0.9982001) q[1];
rz(-pi) q[2];
rz(-2.9473726) q[3];
sx q[3];
rz(-2.9306539) q[3];
sx q[3];
rz(-2.6489182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6151578) q[2];
sx q[2];
rz(-1.7650471) q[2];
sx q[2];
rz(-0.60205013) q[2];
rz(2.708882) q[3];
sx q[3];
rz(-1.5940462) q[3];
sx q[3];
rz(-1.4759147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3848307) q[0];
sx q[0];
rz(-2.6669406) q[0];
sx q[0];
rz(2.5153644) q[0];
rz(-1.2681883) q[1];
sx q[1];
rz(-1.7363997) q[1];
sx q[1];
rz(0.38571206) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6749812) q[0];
sx q[0];
rz(-2.2902596) q[0];
sx q[0];
rz(2.5541259) q[0];
x q[1];
rz(0.63718225) q[2];
sx q[2];
rz(-2.4546625) q[2];
sx q[2];
rz(1.0174583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9358237) q[1];
sx q[1];
rz(-2.9301634) q[1];
sx q[1];
rz(-2.8185259) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3703501) q[3];
sx q[3];
rz(-0.55209898) q[3];
sx q[3];
rz(0.31262661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4422153) q[2];
sx q[2];
rz(-1.135004) q[2];
sx q[2];
rz(-2.8507612) q[2];
rz(-1.95131) q[3];
sx q[3];
rz(-2.5255327) q[3];
sx q[3];
rz(-1.0440089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46665835) q[0];
sx q[0];
rz(-0.062677296) q[0];
sx q[0];
rz(1.689893) q[0];
rz(-0.85707227) q[1];
sx q[1];
rz(-1.8430201) q[1];
sx q[1];
rz(2.7136386) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9687405) q[0];
sx q[0];
rz(-1.4442872) q[0];
sx q[0];
rz(-1.1304024) q[0];
x q[1];
rz(0.54938118) q[2];
sx q[2];
rz(-2.3231988) q[2];
sx q[2];
rz(0.1521509) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9792069) q[1];
sx q[1];
rz(-2.2507994) q[1];
sx q[1];
rz(1.1080145) q[1];
rz(-2.4492743) q[3];
sx q[3];
rz(-2.025617) q[3];
sx q[3];
rz(-2.2934283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9383303) q[2];
sx q[2];
rz(-2.0945956) q[2];
sx q[2];
rz(-0.14673512) q[2];
rz(0.19715582) q[3];
sx q[3];
rz(-1.6517755) q[3];
sx q[3];
rz(2.3728235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41106221) q[0];
sx q[0];
rz(-2.125232) q[0];
sx q[0];
rz(-0.16832571) q[0];
rz(2.7492145) q[1];
sx q[1];
rz(-0.43586755) q[1];
sx q[1];
rz(-1.6900774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0453059) q[0];
sx q[0];
rz(-1.1340965) q[0];
sx q[0];
rz(3.1179881) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8807202) q[2];
sx q[2];
rz(-2.2745273) q[2];
sx q[2];
rz(1.0359302) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5785006) q[1];
sx q[1];
rz(-0.76116409) q[1];
sx q[1];
rz(2.9369686) q[1];
x q[2];
rz(2.6951287) q[3];
sx q[3];
rz(-2.3013448) q[3];
sx q[3];
rz(-1.8346759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3671941) q[2];
sx q[2];
rz(-2.6937679) q[2];
sx q[2];
rz(1.1773342) q[2];
rz(-2.2199953) q[3];
sx q[3];
rz(-2.2955743) q[3];
sx q[3];
rz(-0.53068501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0717936) q[0];
sx q[0];
rz(-1.8853747) q[0];
sx q[0];
rz(-1.0569093) q[0];
rz(0.22843703) q[1];
sx q[1];
rz(-2.3616796) q[1];
sx q[1];
rz(0.25100073) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6901272) q[0];
sx q[0];
rz(-1.8553892) q[0];
sx q[0];
rz(1.837912) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86313451) q[2];
sx q[2];
rz(-0.8388817) q[2];
sx q[2];
rz(1.2108491) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38591641) q[1];
sx q[1];
rz(-1.4078914) q[1];
sx q[1];
rz(-1.4336804) q[1];
rz(-pi) q[2];
rz(-1.1993221) q[3];
sx q[3];
rz(-1.087093) q[3];
sx q[3];
rz(-0.35292816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.935427) q[2];
sx q[2];
rz(-2.3212104) q[2];
sx q[2];
rz(2.7247735) q[2];
rz(-0.56944141) q[3];
sx q[3];
rz(-2.0408401) q[3];
sx q[3];
rz(0.69863629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5702629) q[0];
sx q[0];
rz(-1.8710192) q[0];
sx q[0];
rz(-0.50462333) q[0];
rz(0.2568256) q[1];
sx q[1];
rz(-2.3378614) q[1];
sx q[1];
rz(-1.1508734) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52792943) q[0];
sx q[0];
rz(-1.7759893) q[0];
sx q[0];
rz(2.0733207) q[0];
rz(-0.35607349) q[2];
sx q[2];
rz(-2.1752173) q[2];
sx q[2];
rz(2.6170237) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0432277) q[1];
sx q[1];
rz(-2.4154764) q[1];
sx q[1];
rz(0.67732201) q[1];
rz(-pi) q[2];
rz(-1.7289812) q[3];
sx q[3];
rz(-1.158466) q[3];
sx q[3];
rz(1.7540501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.63200942) q[2];
sx q[2];
rz(-2.7183967) q[2];
sx q[2];
rz(2.7022341) q[2];
rz(-1.1257233) q[3];
sx q[3];
rz(-1.6042387) q[3];
sx q[3];
rz(1.2996947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-0.87089649) q[0];
sx q[0];
rz(-2.3211711) q[0];
sx q[0];
rz(-0.46723715) q[0];
rz(2.1029419) q[1];
sx q[1];
rz(-2.6324582) q[1];
sx q[1];
rz(1.4899303) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3154513) q[0];
sx q[0];
rz(-0.52567654) q[0];
sx q[0];
rz(-2.7868411) q[0];
rz(1.4760255) q[2];
sx q[2];
rz(-1.4247923) q[2];
sx q[2];
rz(0.80751792) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0831775) q[1];
sx q[1];
rz(-1.5442368) q[1];
sx q[1];
rz(2.3258464) q[1];
x q[2];
rz(2.49754) q[3];
sx q[3];
rz(-2.4706512) q[3];
sx q[3];
rz(-1.80598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6963639) q[2];
sx q[2];
rz(-2.7859521) q[2];
sx q[2];
rz(-1.5544372) q[2];
rz(-2.1571531) q[3];
sx q[3];
rz(-2.2319904) q[3];
sx q[3];
rz(1.4136774) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54870355) q[0];
sx q[0];
rz(-0.84243542) q[0];
sx q[0];
rz(0.62826759) q[0];
rz(0.20333044) q[1];
sx q[1];
rz(-2.0275828) q[1];
sx q[1];
rz(-2.5471953) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83594506) q[0];
sx q[0];
rz(-1.3653269) q[0];
sx q[0];
rz(0.97272689) q[0];
x q[1];
rz(0.35697414) q[2];
sx q[2];
rz(-2.5737615) q[2];
sx q[2];
rz(0.50895509) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0600784) q[1];
sx q[1];
rz(-1.8412672) q[1];
sx q[1];
rz(0.31281506) q[1];
rz(-pi) q[2];
rz(-1.6420308) q[3];
sx q[3];
rz(-2.1826577) q[3];
sx q[3];
rz(-2.9737986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9336885) q[2];
sx q[2];
rz(-0.7242569) q[2];
sx q[2];
rz(0.17624632) q[2];
rz(-1.8267953) q[3];
sx q[3];
rz(-2.3254471) q[3];
sx q[3];
rz(2.3740681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15277302) q[0];
sx q[0];
rz(-1.7012699) q[0];
sx q[0];
rz(-1.758601) q[0];
rz(0.11463166) q[1];
sx q[1];
rz(-2.3873867) q[1];
sx q[1];
rz(-0.65232123) q[1];
rz(-2.9822275) q[2];
sx q[2];
rz(-0.2401611) q[2];
sx q[2];
rz(-2.9693749) q[2];
rz(-2.7588941) q[3];
sx q[3];
rz(-1.1094339) q[3];
sx q[3];
rz(-1.4239428) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

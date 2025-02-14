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
rz(-1.3402101) q[0];
sx q[0];
rz(3.4184472) q[0];
sx q[0];
rz(10.463538) q[0];
rz(2.7998595) q[1];
sx q[1];
rz(-0.83655292) q[1];
sx q[1];
rz(-0.41681448) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5688516) q[0];
sx q[0];
rz(-2.4089455) q[0];
sx q[0];
rz(-2.9824663) q[0];
rz(-1.8174581) q[2];
sx q[2];
rz(-1.0283054) q[2];
sx q[2];
rz(-0.88298029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15910251) q[1];
sx q[1];
rz(-1.0034518) q[1];
sx q[1];
rz(-2.2654387) q[1];
rz(-pi) q[2];
rz(2.9800426) q[3];
sx q[3];
rz(-1.6027502) q[3];
sx q[3];
rz(2.3475501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9716399) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(1.2539585) q[2];
rz(-1.5818671) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(1.1221277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9480243) q[0];
sx q[0];
rz(-1.7727611) q[0];
sx q[0];
rz(-2.3728306) q[0];
rz(1.7747152) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(2.1034525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6383253) q[0];
sx q[0];
rz(-1.5621788) q[0];
sx q[0];
rz(1.5805954) q[0];
rz(-0.72121303) q[2];
sx q[2];
rz(-1.8497457) q[2];
sx q[2];
rz(0.21706377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97538725) q[1];
sx q[1];
rz(-2.176099) q[1];
sx q[1];
rz(-1.1206085) q[1];
rz(-2.059405) q[3];
sx q[3];
rz(-1.5540136) q[3];
sx q[3];
rz(0.29573655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64608964) q[2];
sx q[2];
rz(-1.8547408) q[2];
sx q[2];
rz(-2.3213279) q[2];
rz(-0.25343728) q[3];
sx q[3];
rz(-2.7866252) q[3];
sx q[3];
rz(-0.36111116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28061098) q[0];
sx q[0];
rz(-0.51787037) q[0];
sx q[0];
rz(-2.2542727) q[0];
rz(-0.53030983) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(1.1393772) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.087505) q[0];
sx q[0];
rz(-1.9158792) q[0];
sx q[0];
rz(1.4201384) q[0];
rz(-pi) q[1];
rz(-0.33311756) q[2];
sx q[2];
rz(-1.0651759) q[2];
sx q[2];
rz(-0.62084711) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4344221) q[1];
sx q[1];
rz(-1.5853197) q[1];
sx q[1];
rz(2.8959031) q[1];
x q[2];
rz(-1.5185131) q[3];
sx q[3];
rz(-1.0526592) q[3];
sx q[3];
rz(0.62925807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.016971074) q[2];
sx q[2];
rz(-1.6861702) q[2];
sx q[2];
rz(-0.43928453) q[2];
rz(-0.47075054) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(-1.9973756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73862326) q[0];
sx q[0];
rz(-0.93337494) q[0];
sx q[0];
rz(3.0261107) q[0];
rz(-3.0237517) q[1];
sx q[1];
rz(-1.203048) q[1];
sx q[1];
rz(-1.8353362) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71824008) q[0];
sx q[0];
rz(-2.6246492) q[0];
sx q[0];
rz(0.96036185) q[0];
x q[1];
rz(1.9491862) q[2];
sx q[2];
rz(-2.1956964) q[2];
sx q[2];
rz(-1.6922127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0646806) q[1];
sx q[1];
rz(-2.2222328) q[1];
sx q[1];
rz(0.76131911) q[1];
rz(-pi) q[2];
rz(1.4545069) q[3];
sx q[3];
rz(-0.57580417) q[3];
sx q[3];
rz(2.7888014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.304504) q[2];
sx q[2];
rz(-1.3136761) q[2];
sx q[2];
rz(2.9869249) q[2];
rz(1.6020487) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(0.60607564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57347572) q[0];
sx q[0];
rz(-2.0982168) q[0];
sx q[0];
rz(-2.306275) q[0];
rz(2.4256445) q[1];
sx q[1];
rz(-2.7138111) q[1];
sx q[1];
rz(-3.0221525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12163945) q[0];
sx q[0];
rz(-0.058319133) q[0];
sx q[0];
rz(-1.8644237) q[0];
rz(0.44942707) q[2];
sx q[2];
rz(-1.4327421) q[2];
sx q[2];
rz(1.4269478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9000098) q[1];
sx q[1];
rz(-1.1092343) q[1];
sx q[1];
rz(-0.66630967) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64924134) q[3];
sx q[3];
rz(-0.63707817) q[3];
sx q[3];
rz(-2.0811045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9014088) q[2];
sx q[2];
rz(-0.26075026) q[2];
sx q[2];
rz(0.060997941) q[2];
rz(-3.0933464) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(-2.4129996) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.401684) q[0];
sx q[0];
rz(-2.4329199) q[0];
sx q[0];
rz(-2.0106864) q[0];
rz(-2.6629958) q[1];
sx q[1];
rz(-2.5391948) q[1];
sx q[1];
rz(1.4071677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76001747) q[0];
sx q[0];
rz(-1.572346) q[0];
sx q[0];
rz(0.053561915) q[0];
rz(-pi) q[1];
rz(0.84848225) q[2];
sx q[2];
rz(-1.3568078) q[2];
sx q[2];
rz(-3.1297562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0959655) q[1];
sx q[1];
rz(-1.6842323) q[1];
sx q[1];
rz(1.3413702) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9286372) q[3];
sx q[3];
rz(-1.3939438) q[3];
sx q[3];
rz(0.18421728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5799134) q[2];
sx q[2];
rz(-2.734197) q[2];
sx q[2];
rz(-1.178406) q[2];
rz(-1.0493578) q[3];
sx q[3];
rz(-2.6047843) q[3];
sx q[3];
rz(1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2105763) q[0];
sx q[0];
rz(-1.9518305) q[0];
sx q[0];
rz(1.806102) q[0];
rz(1.8203075) q[1];
sx q[1];
rz(-1.0711461) q[1];
sx q[1];
rz(-0.30805045) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0548693) q[0];
sx q[0];
rz(-0.44741524) q[0];
sx q[0];
rz(-2.6347876) q[0];
rz(0.63703434) q[2];
sx q[2];
rz(-2.5278628) q[2];
sx q[2];
rz(-1.6047275) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8845706) q[1];
sx q[1];
rz(-1.0858469) q[1];
sx q[1];
rz(-0.91241769) q[1];
x q[2];
rz(1.6873463) q[3];
sx q[3];
rz(-2.7532059) q[3];
sx q[3];
rz(-1.6167058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41219741) q[2];
sx q[2];
rz(-2.2008379) q[2];
sx q[2];
rz(-3.0103053) q[2];
rz(2.4042551) q[3];
sx q[3];
rz(-1.3056583) q[3];
sx q[3];
rz(2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3103631) q[0];
sx q[0];
rz(-1.0897626) q[0];
sx q[0];
rz(2.6112153) q[0];
rz(-1.4250379) q[1];
sx q[1];
rz(-1.1385463) q[1];
sx q[1];
rz(2.4200965) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5576241) q[0];
sx q[0];
rz(-2.1707105) q[0];
sx q[0];
rz(-2.1198089) q[0];
rz(-0.4690629) q[2];
sx q[2];
rz(-1.4387759) q[2];
sx q[2];
rz(-1.309883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4911982) q[1];
sx q[1];
rz(-1.7868687) q[1];
sx q[1];
rz(-0.19503991) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0305392) q[3];
sx q[3];
rz(-0.78563443) q[3];
sx q[3];
rz(-1.3291886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43712744) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(1.4777769) q[2];
rz(-1.1635228) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(1.6400953) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32996938) q[0];
sx q[0];
rz(-1.7003308) q[0];
sx q[0];
rz(0.60923088) q[0];
rz(1.5513264) q[1];
sx q[1];
rz(-2.8158999) q[1];
sx q[1];
rz(-1.4564266) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0349755) q[0];
sx q[0];
rz(-1.2143232) q[0];
sx q[0];
rz(-0.26654213) q[0];
rz(-0.48522075) q[2];
sx q[2];
rz(-1.6246129) q[2];
sx q[2];
rz(-2.9973928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2472154) q[1];
sx q[1];
rz(-2.651582) q[1];
sx q[1];
rz(-1.9202597) q[1];
rz(-pi) q[2];
rz(2.5129065) q[3];
sx q[3];
rz(-2.1910724) q[3];
sx q[3];
rz(-1.5076021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7234601) q[2];
sx q[2];
rz(-2.2432566) q[2];
sx q[2];
rz(-0.88031236) q[2];
rz(2.9355925) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(-2.6515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708165) q[0];
sx q[0];
rz(-2.4801065) q[0];
sx q[0];
rz(3.1296375) q[0];
rz(-1.6318343) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(2.5451122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2947114) q[0];
sx q[0];
rz(-2.6718326) q[0];
sx q[0];
rz(0.042094783) q[0];
x q[1];
rz(-0.61052236) q[2];
sx q[2];
rz(-0.43402616) q[2];
sx q[2];
rz(2.4422925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4778127) q[1];
sx q[1];
rz(-0.48250178) q[1];
sx q[1];
rz(-1.4003808) q[1];
rz(-pi) q[2];
rz(0.29774547) q[3];
sx q[3];
rz(-1.4199054) q[3];
sx q[3];
rz(2.6835364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.97800469) q[2];
sx q[2];
rz(-0.96077335) q[2];
sx q[2];
rz(-0.19700024) q[2];
rz(-1.9893076) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(-0.44767374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86354179) q[0];
sx q[0];
rz(-2.1320237) q[0];
sx q[0];
rz(-0.19436819) q[0];
rz(2.9327783) q[1];
sx q[1];
rz(-1.5369692) q[1];
sx q[1];
rz(2.1280638) q[1];
rz(1.456719) q[2];
sx q[2];
rz(-2.4488505) q[2];
sx q[2];
rz(1.746576) q[2];
rz(-2.267425) q[3];
sx q[3];
rz(-0.71577358) q[3];
sx q[3];
rz(-1.7795455) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

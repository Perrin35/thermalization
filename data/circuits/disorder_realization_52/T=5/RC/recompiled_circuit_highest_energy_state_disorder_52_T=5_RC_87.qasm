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
rz(1.7088543) q[0];
sx q[0];
rz(6.611293) q[0];
sx q[0];
rz(10.264504) q[0];
rz(-1.7307164) q[1];
sx q[1];
rz(-1.6003992) q[1];
sx q[1];
rz(2.3720001) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1969944) q[0];
sx q[0];
rz(-0.086769484) q[0];
sx q[0];
rz(-1.9831884) q[0];
rz(-pi) q[1];
rz(3.1159471) q[2];
sx q[2];
rz(-0.8249976) q[2];
sx q[2];
rz(0.96482575) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.85629) q[1];
sx q[1];
rz(-1.5714312) q[1];
sx q[1];
rz(-1.2190422) q[1];
rz(-pi) q[2];
rz(-0.64842865) q[3];
sx q[3];
rz(-0.75639137) q[3];
sx q[3];
rz(1.7738916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2644234) q[2];
sx q[2];
rz(-1.2494272) q[2];
sx q[2];
rz(1.9682311) q[2];
rz(2.7133283) q[3];
sx q[3];
rz(-0.56848017) q[3];
sx q[3];
rz(0.65304023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1332755) q[0];
sx q[0];
rz(-2.7011217) q[0];
sx q[0];
rz(1.0622729) q[0];
rz(-1.5509037) q[1];
sx q[1];
rz(-0.56113243) q[1];
sx q[1];
rz(2.0468457) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4832981) q[0];
sx q[0];
rz(-1.6219369) q[0];
sx q[0];
rz(-1.8258894) q[0];
rz(2.6954912) q[2];
sx q[2];
rz(-0.94233905) q[2];
sx q[2];
rz(2.7035509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1909263) q[1];
sx q[1];
rz(-2.6537173) q[1];
sx q[1];
rz(2.4937954) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3571489) q[3];
sx q[3];
rz(-2.3348138) q[3];
sx q[3];
rz(-0.90863228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4636479) q[2];
sx q[2];
rz(-2.8539168) q[2];
sx q[2];
rz(-0.90052432) q[2];
rz(2.1053704) q[3];
sx q[3];
rz(-3.0039054) q[3];
sx q[3];
rz(-0.010995939) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66505945) q[0];
sx q[0];
rz(-2.0158975) q[0];
sx q[0];
rz(1.6735459) q[0];
rz(-2.2131069) q[1];
sx q[1];
rz(-2.1378345) q[1];
sx q[1];
rz(1.4220062) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4741459) q[0];
sx q[0];
rz(-2.8774539) q[0];
sx q[0];
rz(0.35281065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.456383) q[2];
sx q[2];
rz(-2.1985719) q[2];
sx q[2];
rz(0.9762828) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3513438) q[1];
sx q[1];
rz(-2.4741052) q[1];
sx q[1];
rz(-0.60524596) q[1];
x q[2];
rz(0.57392759) q[3];
sx q[3];
rz(-2.6477154) q[3];
sx q[3];
rz(2.9536794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9906087) q[2];
sx q[2];
rz(-1.7093806) q[2];
sx q[2];
rz(-2.3217311) q[2];
rz(3.0153583) q[3];
sx q[3];
rz(-1.1773959) q[3];
sx q[3];
rz(1.646515) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6821297) q[0];
sx q[0];
rz(-1.4012902) q[0];
sx q[0];
rz(1.6023741) q[0];
rz(1.0914717) q[1];
sx q[1];
rz(-2.3075054) q[1];
sx q[1];
rz(1.1702671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226886) q[0];
sx q[0];
rz(-2.6391374) q[0];
sx q[0];
rz(2.2751341) q[0];
rz(-pi) q[1];
rz(-1.865388) q[2];
sx q[2];
rz(-2.3995993) q[2];
sx q[2];
rz(2.2592253) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.023686743) q[1];
sx q[1];
rz(-1.4906018) q[1];
sx q[1];
rz(-3.0515494) q[1];
rz(-0.50962944) q[3];
sx q[3];
rz(-0.42181236) q[3];
sx q[3];
rz(-0.065936655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.22152659) q[2];
sx q[2];
rz(-0.92851323) q[2];
sx q[2];
rz(-0.083219223) q[2];
rz(-0.65256882) q[3];
sx q[3];
rz(-0.021952732) q[3];
sx q[3];
rz(-2.2799802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49587747) q[0];
sx q[0];
rz(-2.8577514) q[0];
sx q[0];
rz(0.19677095) q[0];
rz(-1.981005) q[1];
sx q[1];
rz(-2.6996758) q[1];
sx q[1];
rz(1.7534076) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8179479) q[0];
sx q[0];
rz(-2.9716688) q[0];
sx q[0];
rz(-3.1325045) q[0];
rz(-pi) q[1];
rz(0.20241995) q[2];
sx q[2];
rz(-1.7192063) q[2];
sx q[2];
rz(2.9925516) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39972389) q[1];
sx q[1];
rz(-2.464964) q[1];
sx q[1];
rz(-1.8477738) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9742244) q[3];
sx q[3];
rz(-0.77047548) q[3];
sx q[3];
rz(1.4899173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5768726) q[2];
sx q[2];
rz(-1.2083961) q[2];
sx q[2];
rz(1.2811071) q[2];
rz(-0.34230226) q[3];
sx q[3];
rz(-0.050693158) q[3];
sx q[3];
rz(0.92882338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7814735) q[0];
sx q[0];
rz(-2.0442648) q[0];
sx q[0];
rz(-2.1088364) q[0];
rz(-0.67689854) q[1];
sx q[1];
rz(-0.38336661) q[1];
sx q[1];
rz(-2.3597609) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9932258) q[0];
sx q[0];
rz(-1.9253823) q[0];
sx q[0];
rz(-0.605159) q[0];
rz(-pi) q[1];
x q[1];
rz(1.621945) q[2];
sx q[2];
rz(-1.5384073) q[2];
sx q[2];
rz(0.50078228) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9986753) q[1];
sx q[1];
rz(-1.8774722) q[1];
sx q[1];
rz(0.59345133) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12198066) q[3];
sx q[3];
rz(-2.8013419) q[3];
sx q[3];
rz(1.550479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3013762) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(-2.2678579) q[2];
rz(-0.78911632) q[3];
sx q[3];
rz(-1.5566166) q[3];
sx q[3];
rz(-2.6553787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18272045) q[0];
sx q[0];
rz(-0.046165753) q[0];
sx q[0];
rz(-0.090855457) q[0];
rz(-2.5317522) q[1];
sx q[1];
rz(-1.5515386) q[1];
sx q[1];
rz(-0.59744936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081735805) q[0];
sx q[0];
rz(-1.4422405) q[0];
sx q[0];
rz(-0.14692326) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2944974) q[2];
sx q[2];
rz(-2.6604746) q[2];
sx q[2];
rz(1.4916949) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8272463) q[1];
sx q[1];
rz(-1.1097849) q[1];
sx q[1];
rz(-2.6516312) q[1];
x q[2];
rz(-1.6555637) q[3];
sx q[3];
rz(-2.7206507) q[3];
sx q[3];
rz(1.3011025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37683836) q[2];
sx q[2];
rz(-2.59616) q[2];
sx q[2];
rz(-3.0010014) q[2];
rz(1.2815255) q[3];
sx q[3];
rz(-1.8257273) q[3];
sx q[3];
rz(-0.52711058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.918688) q[0];
sx q[0];
rz(-1.0162901) q[0];
sx q[0];
rz(-1.6131529) q[0];
rz(2.3279922) q[1];
sx q[1];
rz(-3.0366812) q[1];
sx q[1];
rz(2.334107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7225689) q[0];
sx q[0];
rz(-1.3815855) q[0];
sx q[0];
rz(1.7721227) q[0];
rz(-2.0139461) q[2];
sx q[2];
rz(-1.4535744) q[2];
sx q[2];
rz(1.6991311) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85881104) q[1];
sx q[1];
rz(-0.75560299) q[1];
sx q[1];
rz(1.0048466) q[1];
rz(-pi) q[2];
rz(-2.7164338) q[3];
sx q[3];
rz(-1.0711373) q[3];
sx q[3];
rz(2.5806346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7878824) q[2];
sx q[2];
rz(-1.3545802) q[2];
sx q[2];
rz(2.9969969) q[2];
rz(2.0374129) q[3];
sx q[3];
rz(-0.2140597) q[3];
sx q[3];
rz(-2.1000699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5398194) q[0];
sx q[0];
rz(-1.1310534) q[0];
sx q[0];
rz(-2.5850776) q[0];
rz(0.76983184) q[1];
sx q[1];
rz(-0.12431215) q[1];
sx q[1];
rz(-0.46129033) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4249864) q[0];
sx q[0];
rz(-1.5231774) q[0];
sx q[0];
rz(-0.99999962) q[0];
x q[1];
rz(2.2940852) q[2];
sx q[2];
rz(-2.2605148) q[2];
sx q[2];
rz(2.412852) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6016695) q[1];
sx q[1];
rz(-0.8569255) q[1];
sx q[1];
rz(-1.1779835) q[1];
rz(-pi) q[2];
rz(-1.1848524) q[3];
sx q[3];
rz(-2.1084774) q[3];
sx q[3];
rz(-2.4019634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7028659) q[2];
sx q[2];
rz(-2.8361969) q[2];
sx q[2];
rz(0.76496441) q[2];
rz(1.2139828) q[3];
sx q[3];
rz(-0.58911222) q[3];
sx q[3];
rz(-0.37863076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.42534378) q[0];
sx q[0];
rz(-0.049976293) q[0];
sx q[0];
rz(0.36323994) q[0];
rz(-1.4092457) q[1];
sx q[1];
rz(-2.1556518) q[1];
sx q[1];
rz(-2.9149616) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55174819) q[0];
sx q[0];
rz(-1.706011) q[0];
sx q[0];
rz(-1.5677794) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74966263) q[2];
sx q[2];
rz(-1.0790624) q[2];
sx q[2];
rz(-2.3305388) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2001463) q[1];
sx q[1];
rz(-0.90231239) q[1];
sx q[1];
rz(2.723395) q[1];
rz(-pi) q[2];
rz(-1.5814621) q[3];
sx q[3];
rz(-1.9086325) q[3];
sx q[3];
rz(0.84381553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54214415) q[2];
sx q[2];
rz(-0.52079529) q[2];
sx q[2];
rz(0.64860541) q[2];
rz(-3.0829698) q[3];
sx q[3];
rz(-0.62871814) q[3];
sx q[3];
rz(2.4962943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7764353) q[0];
sx q[0];
rz(-1.2189652) q[0];
sx q[0];
rz(-0.49268876) q[0];
rz(-2.0877214) q[1];
sx q[1];
rz(-2.5851879) q[1];
sx q[1];
rz(-2.7256706) q[1];
rz(1.6719978) q[2];
sx q[2];
rz(-0.23739021) q[2];
sx q[2];
rz(-1.8802735) q[2];
rz(-0.875474) q[3];
sx q[3];
rz(-1.3093822) q[3];
sx q[3];
rz(1.1817684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(0.93267814) q[0];
sx q[0];
rz(-2.7490766) q[0];
sx q[0];
rz(2.0445332) q[0];
rz(1.2708083) q[1];
sx q[1];
rz(-1.1975809) q[1];
sx q[1];
rz(1.4374011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2559833) q[0];
sx q[0];
rz(-2.3893111) q[0];
sx q[0];
rz(-0.52623574) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6439866) q[2];
sx q[2];
rz(-1.8022259) q[2];
sx q[2];
rz(-2.481593) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4258051) q[1];
sx q[1];
rz(-2.1658848) q[1];
sx q[1];
rz(-2.7462105) q[1];
rz(0.64793628) q[3];
sx q[3];
rz(-1.5790695) q[3];
sx q[3];
rz(-3.0870761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8561594) q[2];
sx q[2];
rz(-0.23468748) q[2];
sx q[2];
rz(-2.6402546) q[2];
rz(-1.8173789) q[3];
sx q[3];
rz(-0.54558498) q[3];
sx q[3];
rz(-1.6747624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.685574) q[0];
sx q[0];
rz(-0.40559232) q[0];
sx q[0];
rz(1.4612041) q[0];
rz(-0.17924084) q[1];
sx q[1];
rz(-0.76911175) q[1];
sx q[1];
rz(2.5040212) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3831275) q[0];
sx q[0];
rz(-0.33184856) q[0];
sx q[0];
rz(-2.3735574) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.045863002) q[2];
sx q[2];
rz(-1.8872607) q[2];
sx q[2];
rz(-2.4212011) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0936446) q[1];
sx q[1];
rz(-0.94442299) q[1];
sx q[1];
rz(0.44194024) q[1];
rz(-pi) q[2];
rz(2.4781391) q[3];
sx q[3];
rz(-1.2689948) q[3];
sx q[3];
rz(0.38340195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8403988) q[2];
sx q[2];
rz(-2.0638778) q[2];
sx q[2];
rz(-3.0974498) q[2];
rz(0.61331493) q[3];
sx q[3];
rz(-1.3215348) q[3];
sx q[3];
rz(0.80673748) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.089207) q[0];
sx q[0];
rz(-3.0969924) q[0];
sx q[0];
rz(1.3432304) q[0];
rz(-1.7228458) q[1];
sx q[1];
rz(-2.2371465) q[1];
sx q[1];
rz(-2.338063) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65261025) q[0];
sx q[0];
rz(-1.6296415) q[0];
sx q[0];
rz(-0.68025689) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5577689) q[2];
sx q[2];
rz(-1.1808504) q[2];
sx q[2];
rz(-0.7240629) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60119286) q[1];
sx q[1];
rz(-1.4638888) q[1];
sx q[1];
rz(-1.8501297) q[1];
rz(-pi) q[2];
rz(-2.6355631) q[3];
sx q[3];
rz(-1.2057736) q[3];
sx q[3];
rz(-2.6377691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.558305) q[2];
sx q[2];
rz(-2.1896157) q[2];
sx q[2];
rz(-1.0179016) q[2];
rz(1.2244276) q[3];
sx q[3];
rz(-1.3320351) q[3];
sx q[3];
rz(2.0126191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1430436) q[0];
sx q[0];
rz(-0.7989378) q[0];
sx q[0];
rz(2.1813188) q[0];
rz(0.1943365) q[1];
sx q[1];
rz(-1.4413709) q[1];
sx q[1];
rz(-0.0050553102) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.968648) q[0];
sx q[0];
rz(-0.93679777) q[0];
sx q[0];
rz(0.89374884) q[0];
rz(-1.5326535) q[2];
sx q[2];
rz(-1.3873709) q[2];
sx q[2];
rz(-2.9596552) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21987629) q[1];
sx q[1];
rz(-1.7895234) q[1];
sx q[1];
rz(1.686386) q[1];
rz(-pi) q[2];
rz(-2.9937135) q[3];
sx q[3];
rz(-1.0564466) q[3];
sx q[3];
rz(-1.0837931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2527577) q[2];
sx q[2];
rz(-1.5965261) q[2];
sx q[2];
rz(1.6952391) q[2];
rz(-1.3085922) q[3];
sx q[3];
rz(-1.3819709) q[3];
sx q[3];
rz(0.60477177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5204891) q[0];
sx q[0];
rz(-2.4172754) q[0];
sx q[0];
rz(-2.5526168) q[0];
rz(0.0079872459) q[1];
sx q[1];
rz(-2.0365448) q[1];
sx q[1];
rz(1.084682) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5552591) q[0];
sx q[0];
rz(-1.8249092) q[0];
sx q[0];
rz(1.8122559) q[0];
x q[1];
rz(-0.029726278) q[2];
sx q[2];
rz(-2.8619302) q[2];
sx q[2];
rz(-0.047911876) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1639203) q[1];
sx q[1];
rz(-2.2798988) q[1];
sx q[1];
rz(-0.52859882) q[1];
rz(-pi) q[2];
rz(0.9081697) q[3];
sx q[3];
rz(-1.623933) q[3];
sx q[3];
rz(-2.4897864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7320554) q[2];
sx q[2];
rz(-2.1246702) q[2];
sx q[2];
rz(-1.3746877) q[2];
rz(3.0459259) q[3];
sx q[3];
rz(-1.5868264) q[3];
sx q[3];
rz(-2.5110551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.10504) q[0];
sx q[0];
rz(-0.51346546) q[0];
sx q[0];
rz(1.8010944) q[0];
rz(-2.4758677) q[1];
sx q[1];
rz(-1.7981073) q[1];
sx q[1];
rz(-0.7241157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.501717) q[0];
sx q[0];
rz(-2.6611009) q[0];
sx q[0];
rz(2.2697422) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3409136) q[2];
sx q[2];
rz(-0.73298798) q[2];
sx q[2];
rz(1.6782325) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16038475) q[1];
sx q[1];
rz(-1.1120468) q[1];
sx q[1];
rz(-1.4190516) q[1];
rz(-1.4744711) q[3];
sx q[3];
rz(-1.5755638) q[3];
sx q[3];
rz(-1.7624439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93806997) q[2];
sx q[2];
rz(-1.6911493) q[2];
sx q[2];
rz(3.0435437) q[2];
rz(-0.34052643) q[3];
sx q[3];
rz(-2.2398658) q[3];
sx q[3];
rz(0.19671973) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53511867) q[0];
sx q[0];
rz(-2.4381194) q[0];
sx q[0];
rz(-0.055543609) q[0];
rz(2.7659888) q[1];
sx q[1];
rz(-0.77508488) q[1];
sx q[1];
rz(-1.6966049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7168769) q[0];
sx q[0];
rz(-2.3622241) q[0];
sx q[0];
rz(-2.9911441) q[0];
rz(-pi) q[1];
rz(-1.7710502) q[2];
sx q[2];
rz(-2.8828388) q[2];
sx q[2];
rz(0.090752964) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6003528) q[1];
sx q[1];
rz(-1.9484047) q[1];
sx q[1];
rz(1.5560335) q[1];
rz(0.7020601) q[3];
sx q[3];
rz(-2.5581048) q[3];
sx q[3];
rz(-1.7398011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0982509) q[2];
sx q[2];
rz(-1.8378374) q[2];
sx q[2];
rz(-1.8275758) q[2];
rz(-0.029953778) q[3];
sx q[3];
rz(-1.6311389) q[3];
sx q[3];
rz(-2.8225074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444262) q[0];
sx q[0];
rz(-0.63969669) q[0];
sx q[0];
rz(1.8901012) q[0];
rz(-0.80998069) q[1];
sx q[1];
rz(-2.4046343) q[1];
sx q[1];
rz(-1.4553778) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7969709) q[0];
sx q[0];
rz(-0.83528548) q[0];
sx q[0];
rz(-1.263144) q[0];
rz(-0.095429777) q[2];
sx q[2];
rz(-0.52431267) q[2];
sx q[2];
rz(1.4316259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4659459) q[1];
sx q[1];
rz(-2.1928804) q[1];
sx q[1];
rz(0.81333604) q[1];
rz(-2.3076644) q[3];
sx q[3];
rz(-1.0627316) q[3];
sx q[3];
rz(0.64475393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2927085) q[2];
sx q[2];
rz(-1.3564738) q[2];
sx q[2];
rz(-1.1513618) q[2];
rz(-0.014160841) q[3];
sx q[3];
rz(-1.7834981) q[3];
sx q[3];
rz(-1.8748698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8779811) q[0];
sx q[0];
rz(-0.71313715) q[0];
sx q[0];
rz(-0.94733316) q[0];
rz(2.5454648) q[1];
sx q[1];
rz(-1.4214186) q[1];
sx q[1];
rz(-1.5790342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9536586) q[0];
sx q[0];
rz(-0.74043027) q[0];
sx q[0];
rz(-2.8859167) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30513341) q[2];
sx q[2];
rz(-0.49499475) q[2];
sx q[2];
rz(-2.309805) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71982518) q[1];
sx q[1];
rz(-0.9523069) q[1];
sx q[1];
rz(2.5773125) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4189255) q[3];
sx q[3];
rz(-2.8407466) q[3];
sx q[3];
rz(2.2900555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8922213) q[2];
sx q[2];
rz(-2.5669079) q[2];
sx q[2];
rz(0.48386827) q[2];
rz(0.59746915) q[3];
sx q[3];
rz(-1.6941083) q[3];
sx q[3];
rz(-0.57815236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4018965) q[0];
sx q[0];
rz(-1.551349) q[0];
sx q[0];
rz(0.36648146) q[0];
rz(1.3260427) q[1];
sx q[1];
rz(-0.96249023) q[1];
sx q[1];
rz(2.0749626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90049998) q[0];
sx q[0];
rz(-2.0632732) q[0];
sx q[0];
rz(2.2482613) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5353955) q[2];
sx q[2];
rz(-1.5328888) q[2];
sx q[2];
rz(-0.67680046) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45056709) q[1];
sx q[1];
rz(-0.98916173) q[1];
sx q[1];
rz(-0.29009351) q[1];
rz(0.033081369) q[3];
sx q[3];
rz(-2.0478147) q[3];
sx q[3];
rz(-1.5257924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8804973) q[2];
sx q[2];
rz(-0.77793056) q[2];
sx q[2];
rz(-2.0757389) q[2];
rz(-2.8737658) q[3];
sx q[3];
rz(-1.4466176) q[3];
sx q[3];
rz(-2.6563787) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5601226) q[0];
sx q[0];
rz(-2.1385834) q[0];
sx q[0];
rz(-2.9271097) q[0];
rz(2.8196234) q[1];
sx q[1];
rz(-1.8245158) q[1];
sx q[1];
rz(1.239924) q[1];
rz(-1.5787081) q[2];
sx q[2];
rz(-1.2894165) q[2];
sx q[2];
rz(-3.0647562) q[2];
rz(-2.3278676) q[3];
sx q[3];
rz(-0.56725262) q[3];
sx q[3];
rz(-0.57813416) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

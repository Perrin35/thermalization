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
rz(-1.7041915) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8575246) q[0];
sx q[0];
rz(-1.2204613) q[0];
sx q[0];
rz(0.68035462) q[0];
x q[1];
rz(-0.3008879) q[2];
sx q[2];
rz(-2.8990633) q[2];
sx q[2];
rz(-2.1721942) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7157876) q[1];
sx q[1];
rz(-0.97570786) q[1];
sx q[1];
rz(-0.39538212) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.013707073) q[3];
sx q[3];
rz(-2.4936112) q[3];
sx q[3];
rz(-1.614384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8561594) q[2];
sx q[2];
rz(-0.23468748) q[2];
sx q[2];
rz(-2.6402546) q[2];
rz(1.8173789) q[3];
sx q[3];
rz(-0.54558498) q[3];
sx q[3];
rz(-1.4668303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45601869) q[0];
sx q[0];
rz(-0.40559232) q[0];
sx q[0];
rz(-1.4612041) q[0];
rz(2.9623518) q[1];
sx q[1];
rz(-2.3724809) q[1];
sx q[1];
rz(0.63757149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9623827) q[0];
sx q[0];
rz(-1.80733) q[0];
sx q[0];
rz(1.80577) q[0];
rz(-pi) q[1];
rz(1.4316971) q[2];
sx q[2];
rz(-0.31965986) q[2];
sx q[2];
rz(-2.5676198) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.04794807) q[1];
sx q[1];
rz(-2.1971697) q[1];
sx q[1];
rz(2.6996524) q[1];
rz(-pi) q[2];
rz(-0.66345352) q[3];
sx q[3];
rz(-1.8725978) q[3];
sx q[3];
rz(-0.38340195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8403988) q[2];
sx q[2];
rz(-1.0777148) q[2];
sx q[2];
rz(-3.0974498) q[2];
rz(-0.61331493) q[3];
sx q[3];
rz(-1.8200579) q[3];
sx q[3];
rz(-2.3348552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052385656) q[0];
sx q[0];
rz(-3.0969924) q[0];
sx q[0];
rz(-1.3432304) q[0];
rz(1.4187468) q[1];
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
rz(-0.87063861) q[0];
sx q[0];
rz(-0.89193908) q[0];
sx q[0];
rz(-1.4951597) q[0];
rz(-pi) q[1];
rz(2.0284925) q[2];
sx q[2];
rz(-2.105793) q[2];
sx q[2];
rz(2.0488103) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5281339) q[1];
sx q[1];
rz(-2.8430004) q[1];
sx q[1];
rz(1.9419844) q[1];
rz(-pi) q[2];
rz(1.9827051) q[3];
sx q[3];
rz(-2.0406699) q[3];
sx q[3];
rz(-1.2622693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.558305) q[2];
sx q[2];
rz(-0.95197695) q[2];
sx q[2];
rz(2.1236911) q[2];
rz(-1.917165) q[3];
sx q[3];
rz(-1.8095576) q[3];
sx q[3];
rz(-2.0126191) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99854904) q[0];
sx q[0];
rz(-0.7989378) q[0];
sx q[0];
rz(2.1813188) q[0];
rz(0.1943365) q[1];
sx q[1];
rz(-1.7002218) q[1];
sx q[1];
rz(-3.1365373) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1729447) q[0];
sx q[0];
rz(-0.93679777) q[0];
sx q[0];
rz(-0.89374884) q[0];
x q[1];
rz(0.18355592) q[2];
sx q[2];
rz(-1.5332937) q[2];
sx q[2];
rz(1.3958193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27144602) q[1];
sx q[1];
rz(-0.24695858) q[1];
sx q[1];
rz(-0.4786164) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8258589) q[3];
sx q[3];
rz(-2.6082468) q[3];
sx q[3];
rz(-1.3778138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.888835) q[2];
sx q[2];
rz(-1.5965261) q[2];
sx q[2];
rz(-1.4463536) q[2];
rz(1.8330005) q[3];
sx q[3];
rz(-1.3819709) q[3];
sx q[3];
rz(-2.5368209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5204891) q[0];
sx q[0];
rz(-2.4172754) q[0];
sx q[0];
rz(2.5526168) q[0];
rz(0.0079872459) q[1];
sx q[1];
rz(-1.1050478) q[1];
sx q[1];
rz(2.0569107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046291489) q[0];
sx q[0];
rz(-1.8043552) q[0];
sx q[0];
rz(-0.26136847) q[0];
rz(-3.1118664) q[2];
sx q[2];
rz(-0.27966248) q[2];
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
rz(2.6129938) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.067361319) q[3];
sx q[3];
rz(-2.2323221) q[3];
sx q[3];
rz(-2.1811821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40953723) q[2];
sx q[2];
rz(-2.1246702) q[2];
sx q[2];
rz(-1.3746877) q[2];
rz(0.095666766) q[3];
sx q[3];
rz(-1.5868264) q[3];
sx q[3];
rz(2.5110551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0365527) q[0];
sx q[0];
rz(-2.6281272) q[0];
sx q[0];
rz(1.3404982) q[0];
rz(-0.66572491) q[1];
sx q[1];
rz(-1.7981073) q[1];
sx q[1];
rz(-2.417477) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.501717) q[0];
sx q[0];
rz(-0.48049179) q[0];
sx q[0];
rz(0.87185045) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3409136) q[2];
sx q[2];
rz(-2.4086047) q[2];
sx q[2];
rz(1.6782325) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49289068) q[1];
sx q[1];
rz(-2.6600983) q[1];
sx q[1];
rz(2.8445811) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0047896623) q[3];
sx q[3];
rz(-1.4744722) q[3];
sx q[3];
rz(-0.19118689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2035227) q[2];
sx q[2];
rz(-1.4504434) q[2];
sx q[2];
rz(-0.098048992) q[2];
rz(-0.34052643) q[3];
sx q[3];
rz(-2.2398658) q[3];
sx q[3];
rz(-2.9448729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53511867) q[0];
sx q[0];
rz(-2.4381194) q[0];
sx q[0];
rz(0.055543609) q[0];
rz(-2.7659888) q[1];
sx q[1];
rz(-2.3665078) q[1];
sx q[1];
rz(-1.6966049) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21475131) q[0];
sx q[0];
rz(-2.3390649) q[0];
sx q[0];
rz(-1.4237798) q[0];
rz(-pi) q[1];
rz(1.7710502) q[2];
sx q[2];
rz(-2.8828388) q[2];
sx q[2];
rz(-0.090752964) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58126175) q[1];
sx q[1];
rz(-0.37788299) q[1];
sx q[1];
rz(0.037200971) q[1];
x q[2];
rz(1.1678004) q[3];
sx q[3];
rz(-2.0049565) q[3];
sx q[3];
rz(2.1938965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0982509) q[2];
sx q[2];
rz(-1.3037553) q[2];
sx q[2];
rz(1.8275758) q[2];
rz(-0.029953778) q[3];
sx q[3];
rz(-1.5104537) q[3];
sx q[3];
rz(-0.31908527) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-0.73695838) q[1];
sx q[1];
rz(1.4553778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90235119) q[0];
sx q[0];
rz(-2.3555814) q[0];
sx q[0];
rz(2.8186574) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6258442) q[2];
sx q[2];
rz(-1.0491129) q[2];
sx q[2];
rz(-1.59984) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4474214) q[1];
sx q[1];
rz(-2.2024367) q[1];
sx q[1];
rz(-2.3775435) q[1];
rz(-pi) q[2];
rz(2.2627566) q[3];
sx q[3];
rz(-0.86711001) q[3];
sx q[3];
rz(1.4181942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2927085) q[2];
sx q[2];
rz(-1.7851189) q[2];
sx q[2];
rz(-1.1513618) q[2];
rz(0.014160841) q[3];
sx q[3];
rz(-1.7834981) q[3];
sx q[3];
rz(1.8748698) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8779811) q[0];
sx q[0];
rz(-0.71313715) q[0];
sx q[0];
rz(0.94733316) q[0];
rz(-2.5454648) q[1];
sx q[1];
rz(-1.4214186) q[1];
sx q[1];
rz(1.5790342) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84760188) q[0];
sx q[0];
rz(-2.2819073) q[0];
sx q[0];
rz(1.7979277) q[0];
rz(1.7315699) q[2];
sx q[2];
rz(-2.0410182) q[2];
sx q[2];
rz(1.9660814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49925466) q[1];
sx q[1];
rz(-2.0216536) q[1];
sx q[1];
rz(-0.87085215) q[1];
rz(-pi) q[2];
rz(-0.72266717) q[3];
sx q[3];
rz(-0.30084601) q[3];
sx q[3];
rz(-0.85153714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24937135) q[2];
sx q[2];
rz(-0.5746848) q[2];
sx q[2];
rz(-0.48386827) q[2];
rz(0.59746915) q[3];
sx q[3];
rz(-1.6941083) q[3];
sx q[3];
rz(-0.57815236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4018965) q[0];
sx q[0];
rz(-1.5902436) q[0];
sx q[0];
rz(-0.36648146) q[0];
rz(1.81555) q[1];
sx q[1];
rz(-2.1791024) q[1];
sx q[1];
rz(-1.0666301) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.033762) q[0];
sx q[0];
rz(-0.98563802) q[0];
sx q[0];
rz(2.5385227) q[0];
rz(-pi) q[1];
rz(1.5246806) q[2];
sx q[2];
rz(-0.96509714) q[2];
sx q[2];
rz(2.2213288) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2827833) q[1];
sx q[1];
rz(-1.3294535) q[1];
sx q[1];
rz(2.1721852) q[1];
rz(-pi) q[2];
rz(2.048038) q[3];
sx q[3];
rz(-1.541409) q[3];
sx q[3];
rz(0.060197006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2610953) q[2];
sx q[2];
rz(-2.3636621) q[2];
sx q[2];
rz(1.0658537) q[2];
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
rz(-pi) q[2];
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
rz(0.32196925) q[1];
sx q[1];
rz(-1.3170769) q[1];
sx q[1];
rz(-1.9016686) q[1];
rz(1.5787081) q[2];
sx q[2];
rz(-1.8521761) q[2];
sx q[2];
rz(0.076836486) q[2];
rz(2.0044708) q[3];
sx q[3];
rz(-1.9487582) q[3];
sx q[3];
rz(-1.4760803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-2.5118339) q[0];
sx q[0];
rz(-0.30183733) q[0];
sx q[0];
rz(1.3071741) q[0];
rz(-3.6699927) q[1];
sx q[1];
rz(3.9730605) q[1];
sx q[1];
rz(12.107036) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0889657) q[0];
sx q[0];
rz(-1.3148493) q[0];
sx q[0];
rz(0.050764485) q[0];
x q[1];
rz(2.1856543) q[2];
sx q[2];
rz(-2.8070076) q[2];
sx q[2];
rz(-1.057404) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85831538) q[1];
sx q[1];
rz(-2.1766571) q[1];
sx q[1];
rz(-1.359904) q[1];
rz(0.022597952) q[3];
sx q[3];
rz(-0.90211464) q[3];
sx q[3];
rz(-0.99136664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1186195) q[2];
sx q[2];
rz(-1.8034673) q[2];
sx q[2];
rz(2.5845134) q[2];
rz(-2.6856375) q[3];
sx q[3];
rz(-2.3796701) q[3];
sx q[3];
rz(0.60321155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93422455) q[0];
sx q[0];
rz(-2.4226483) q[0];
sx q[0];
rz(-1.3907322) q[0];
rz(-1.8234183) q[1];
sx q[1];
rz(-1.428182) q[1];
sx q[1];
rz(2.9255829) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0829613) q[0];
sx q[0];
rz(-2.4602232) q[0];
sx q[0];
rz(-2.0860703) q[0];
rz(-pi) q[1];
rz(-2.7223396) q[2];
sx q[2];
rz(-1.6800095) q[2];
sx q[2];
rz(-2.2855503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60571721) q[1];
sx q[1];
rz(-0.89369666) q[1];
sx q[1];
rz(0.89800055) q[1];
rz(-2.1116637) q[3];
sx q[3];
rz(-1.3824302) q[3];
sx q[3];
rz(2.4595367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4231437) q[2];
sx q[2];
rz(-0.37280145) q[2];
sx q[2];
rz(-3.0917523) q[2];
rz(0.54346624) q[3];
sx q[3];
rz(-2.016957) q[3];
sx q[3];
rz(1.0928833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8922888) q[0];
sx q[0];
rz(-1.1019305) q[0];
sx q[0];
rz(0.9683384) q[0];
rz(0.020542055) q[1];
sx q[1];
rz(-1.3803218) q[1];
sx q[1];
rz(-1.7113908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0413301) q[0];
sx q[0];
rz(-1.5484637) q[0];
sx q[0];
rz(2.5002527) q[0];
rz(0.56349748) q[2];
sx q[2];
rz(-1.2545956) q[2];
sx q[2];
rz(1.344156) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8404482) q[1];
sx q[1];
rz(-1.9663875) q[1];
sx q[1];
rz(-1.3592833) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7647721) q[3];
sx q[3];
rz(-1.6018036) q[3];
sx q[3];
rz(-0.12801192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4358383) q[2];
sx q[2];
rz(-0.64283723) q[2];
sx q[2];
rz(1.359681) q[2];
rz(0.5082353) q[3];
sx q[3];
rz(-1.5225007) q[3];
sx q[3];
rz(-1.9967509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.085676) q[0];
sx q[0];
rz(-2.8450232) q[0];
sx q[0];
rz(-2.1348409) q[0];
rz(2.309917) q[1];
sx q[1];
rz(-0.74200231) q[1];
sx q[1];
rz(2.0337909) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2058973) q[0];
sx q[0];
rz(-2.1855121) q[0];
sx q[0];
rz(0.45386916) q[0];
x q[1];
rz(-0.37563373) q[2];
sx q[2];
rz(-1.3357329) q[2];
sx q[2];
rz(-1.7441074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7101871) q[1];
sx q[1];
rz(-1.5905989) q[1];
sx q[1];
rz(-1.2638249) q[1];
rz(1.8516225) q[3];
sx q[3];
rz(-1.1410332) q[3];
sx q[3];
rz(0.44441477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6882249) q[2];
sx q[2];
rz(-1.0658762) q[2];
sx q[2];
rz(-0.31309703) q[2];
rz(-2.744216) q[3];
sx q[3];
rz(-1.671096) q[3];
sx q[3];
rz(-2.9562922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6322286) q[0];
sx q[0];
rz(-0.86357421) q[0];
sx q[0];
rz(-2.2160227) q[0];
rz(-0.46982345) q[1];
sx q[1];
rz(-1.8269822) q[1];
sx q[1];
rz(-1.0161317) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39611577) q[0];
sx q[0];
rz(-1.1594311) q[0];
sx q[0];
rz(1.0883254) q[0];
rz(1.2193331) q[2];
sx q[2];
rz(-1.8242852) q[2];
sx q[2];
rz(-0.070310449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4279856) q[1];
sx q[1];
rz(-1.813289) q[1];
sx q[1];
rz(0.20771435) q[1];
rz(1.8941325) q[3];
sx q[3];
rz(-1.5405021) q[3];
sx q[3];
rz(-0.39112511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4887345) q[2];
sx q[2];
rz(-0.32565871) q[2];
sx q[2];
rz(2.6540836) q[2];
rz(2.2440535) q[3];
sx q[3];
rz(-1.8833912) q[3];
sx q[3];
rz(-2.0770309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2010736) q[0];
sx q[0];
rz(-2.2036393) q[0];
sx q[0];
rz(-2.4679389) q[0];
rz(-1.2334088) q[1];
sx q[1];
rz(-2.1234832) q[1];
sx q[1];
rz(2.3755551) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29391321) q[0];
sx q[0];
rz(-1.3013757) q[0];
sx q[0];
rz(-2.1015757) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38389194) q[2];
sx q[2];
rz(-2.9523179) q[2];
sx q[2];
rz(1.6757019) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1704857) q[1];
sx q[1];
rz(-2.0506713) q[1];
sx q[1];
rz(2.9779153) q[1];
rz(-pi) q[2];
rz(0.30575846) q[3];
sx q[3];
rz(-1.7019203) q[3];
sx q[3];
rz(-2.6763889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.172714) q[2];
sx q[2];
rz(-1.6438899) q[2];
sx q[2];
rz(-0.76095757) q[2];
rz(2.739665) q[3];
sx q[3];
rz(-2.8765078) q[3];
sx q[3];
rz(0.5886122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6589979) q[0];
sx q[0];
rz(-0.75356475) q[0];
sx q[0];
rz(-1.6931417) q[0];
rz(-0.50085577) q[1];
sx q[1];
rz(-1.294699) q[1];
sx q[1];
rz(-0.81333152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017352176) q[0];
sx q[0];
rz(-0.4158786) q[0];
sx q[0];
rz(-2.7680567) q[0];
x q[1];
rz(-3.0221239) q[2];
sx q[2];
rz(-2.1158233) q[2];
sx q[2];
rz(1.8785005) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8706559) q[1];
sx q[1];
rz(-1.7084048) q[1];
sx q[1];
rz(0.44160053) q[1];
x q[2];
rz(2.0175948) q[3];
sx q[3];
rz(-0.13544336) q[3];
sx q[3];
rz(-0.14458421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27998754) q[2];
sx q[2];
rz(-0.6051175) q[2];
sx q[2];
rz(-0.32312265) q[2];
rz(1.7396287) q[3];
sx q[3];
rz(-1.8596545) q[3];
sx q[3];
rz(-1.3824979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0793656) q[0];
sx q[0];
rz(-2.0923738) q[0];
sx q[0];
rz(0.76612377) q[0];
rz(-0.8194204) q[1];
sx q[1];
rz(-1.0304281) q[1];
sx q[1];
rz(1.5310418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19056828) q[0];
sx q[0];
rz(-1.5026717) q[0];
sx q[0];
rz(0.22916746) q[0];
x q[1];
rz(0.65532897) q[2];
sx q[2];
rz(-1.3957784) q[2];
sx q[2];
rz(-1.9536215) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30554015) q[1];
sx q[1];
rz(-2.17872) q[1];
sx q[1];
rz(-2.4120283) q[1];
rz(-2.2977912) q[3];
sx q[3];
rz(-0.42359023) q[3];
sx q[3];
rz(-2.9543878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3766342) q[2];
sx q[2];
rz(-2.9464293) q[2];
sx q[2];
rz(-0.39806077) q[2];
rz(2.5352488) q[3];
sx q[3];
rz(-1.5747993) q[3];
sx q[3];
rz(-2.2280367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.23676087) q[0];
sx q[0];
rz(-2.8396711) q[0];
sx q[0];
rz(2.6171369) q[0];
rz(-0.66871387) q[1];
sx q[1];
rz(-1.6122183) q[1];
sx q[1];
rz(-1.761577) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4008025) q[0];
sx q[0];
rz(-2.5668416) q[0];
sx q[0];
rz(1.4905592) q[0];
x q[1];
rz(-0.30461208) q[2];
sx q[2];
rz(-1.1122983) q[2];
sx q[2];
rz(0.079041399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.804033) q[1];
sx q[1];
rz(-1.9849249) q[1];
sx q[1];
rz(2.0284589) q[1];
rz(-0.59137592) q[3];
sx q[3];
rz(-0.61695652) q[3];
sx q[3];
rz(-0.57720952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1606719) q[2];
sx q[2];
rz(-0.74690861) q[2];
sx q[2];
rz(-0.46372947) q[2];
rz(1.4240228) q[3];
sx q[3];
rz(-1.0878891) q[3];
sx q[3];
rz(3.1080642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5793107) q[0];
sx q[0];
rz(-0.71826851) q[0];
sx q[0];
rz(-0.41807362) q[0];
rz(-2.383291) q[1];
sx q[1];
rz(-0.98379358) q[1];
sx q[1];
rz(-1.2850579) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67659527) q[0];
sx q[0];
rz(-0.9228068) q[0];
sx q[0];
rz(1.4761488) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1056771) q[2];
sx q[2];
rz(-0.67567247) q[2];
sx q[2];
rz(3.1356406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5922136) q[1];
sx q[1];
rz(-2.5587808) q[1];
sx q[1];
rz(0.35415502) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.140652) q[3];
sx q[3];
rz(-2.7469198) q[3];
sx q[3];
rz(-2.3230524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9717504) q[2];
sx q[2];
rz(-1.0047793) q[2];
sx q[2];
rz(2.8954835) q[2];
rz(-1.8698112) q[3];
sx q[3];
rz(-1.5747986) q[3];
sx q[3];
rz(-1.7790214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9825738) q[0];
sx q[0];
rz(-1.4763426) q[0];
sx q[0];
rz(2.939298) q[0];
rz(-2.5166439) q[1];
sx q[1];
rz(-0.85660558) q[1];
sx q[1];
rz(0.4013335) q[1];
rz(-0.18904674) q[2];
sx q[2];
rz(-1.3658267) q[2];
sx q[2];
rz(1.7455802) q[2];
rz(2.940964) q[3];
sx q[3];
rz(-1.8798141) q[3];
sx q[3];
rz(0.78165913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

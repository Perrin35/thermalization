OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971561) q[0];
sx q[0];
rz(-1.7261788) q[0];
sx q[0];
rz(-0.42874254) q[0];
x q[1];
rz(0.46967311) q[2];
sx q[2];
rz(-2.854752) q[2];
sx q[2];
rz(1.6490205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2130148) q[1];
sx q[1];
rz(-1.9933812) q[1];
sx q[1];
rz(-2.1133452) q[1];
rz(0.10098884) q[3];
sx q[3];
rz(-1.0102934) q[3];
sx q[3];
rz(-1.8141754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1564864) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1433379) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(0.026219333) q[0];
rz(-1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(0.96347934) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.11461) q[0];
sx q[0];
rz(-2.5275505) q[0];
sx q[0];
rz(-1.5747889) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0019148) q[2];
sx q[2];
rz(-1.2895673) q[2];
sx q[2];
rz(2.1525454) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7746437) q[1];
sx q[1];
rz(-2.0683214) q[1];
sx q[1];
rz(-0.36555396) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2198592) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(-1.9333145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6271237) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(2.3441558) q[0];
rz(1.047661) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(2.581596) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0021734) q[0];
sx q[0];
rz(-1.9529441) q[0];
sx q[0];
rz(0.21811534) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5494924) q[2];
sx q[2];
rz(-1.8739803) q[2];
sx q[2];
rz(1.6456749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8078976) q[1];
sx q[1];
rz(-1.9296608) q[1];
sx q[1];
rz(-1.7829423) q[1];
rz(1.0873763) q[3];
sx q[3];
rz(-1.1421575) q[3];
sx q[3];
rz(-2.7408858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75227633) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(0.17253549) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.003222) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(1.2987312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7553058) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(-1.132071) q[0];
rz(0.0012871731) q[2];
sx q[2];
rz(-2.3372071) q[2];
sx q[2];
rz(0.019891642) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5609834) q[1];
sx q[1];
rz(-1.0099851) q[1];
sx q[1];
rz(-0.22052712) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8989765) q[3];
sx q[3];
rz(-1.4644074) q[3];
sx q[3];
rz(1.002841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(0.38468012) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(-1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5383179) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(-1.3866562) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(-0.2968266) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291527) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(-2.5705283) q[0];
rz(-pi) q[1];
rz(-0.16634059) q[2];
sx q[2];
rz(-2.3918249) q[2];
sx q[2];
rz(-1.9094085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0592812) q[1];
sx q[1];
rz(-1.1059522) q[1];
sx q[1];
rz(-2.6962198) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4279168) q[3];
sx q[3];
rz(-1.5034961) q[3];
sx q[3];
rz(-1.0825368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7632873) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(0.39247593) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0157938) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(0.75138599) q[0];
rz(-1.3279351) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(-2.5352535) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7327001) q[0];
sx q[0];
rz(-3.0928287) q[0];
sx q[0];
rz(2.7932037) q[0];
rz(-pi) q[1];
rz(0.41646429) q[2];
sx q[2];
rz(-2.205924) q[2];
sx q[2];
rz(-0.093402775) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9669173) q[1];
sx q[1];
rz(-1.3904018) q[1];
sx q[1];
rz(1.0134407) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0665928) q[3];
sx q[3];
rz(-1.787775) q[3];
sx q[3];
rz(1.0523588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5027344) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(1.139337) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(-0.50810057) q[0];
rz(-1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(0.79024822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19827851) q[0];
sx q[0];
rz(-0.6753079) q[0];
sx q[0];
rz(0.4839464) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5590645) q[2];
sx q[2];
rz(-1.3571697) q[2];
sx q[2];
rz(-2.8649883) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33752791) q[1];
sx q[1];
rz(-1.84066) q[1];
sx q[1];
rz(-2.7359664) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1931476) q[3];
sx q[3];
rz(-2.1834063) q[3];
sx q[3];
rz(0.8876422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(1.4132168) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4247894) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.75288) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1202639) q[0];
sx q[0];
rz(-1.6006032) q[0];
sx q[0];
rz(-2.1304312) q[0];
rz(2.4829364) q[2];
sx q[2];
rz(-2.5052862) q[2];
sx q[2];
rz(-0.051740019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74275201) q[1];
sx q[1];
rz(-1.5594149) q[1];
sx q[1];
rz(-0.091817261) q[1];
rz(0.29626131) q[3];
sx q[3];
rz(-2.1564266) q[3];
sx q[3];
rz(-1.1405917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1887112) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(0.88225538) q[2];
rz(1.7404209) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27154487) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.4260938) q[0];
rz(-0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-0.55823278) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4687913) q[0];
sx q[0];
rz(-1.1835915) q[0];
sx q[0];
rz(2.0331435) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4633281) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(-1.4449643) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8856814) q[1];
sx q[1];
rz(-2.3951888) q[1];
sx q[1];
rz(-1.6649151) q[1];
rz(-pi) q[2];
rz(-1.3853119) q[3];
sx q[3];
rz(-2.0501325) q[3];
sx q[3];
rz(-2.6244147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2925064) q[2];
sx q[2];
rz(-1.2693274) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.50487173) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(-0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-2.6182981) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4921724) q[0];
sx q[0];
rz(-1.5086552) q[0];
sx q[0];
rz(-0.6092351) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5027572) q[2];
sx q[2];
rz(-2.5685852) q[2];
sx q[2];
rz(0.21493658) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0303505) q[1];
sx q[1];
rz(-0.48644629) q[1];
sx q[1];
rz(0.72279795) q[1];
rz(-pi) q[2];
rz(0.2449805) q[3];
sx q[3];
rz(-1.7953201) q[3];
sx q[3];
rz(-2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(-2.8722897) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(-2.0555029) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8158648) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(3.1303828) q[2];
sx q[2];
rz(-1.3080454) q[2];
sx q[2];
rz(-1.0946454) q[2];
rz(0.23444093) q[3];
sx q[3];
rz(-0.81080484) q[3];
sx q[3];
rz(-0.57639359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

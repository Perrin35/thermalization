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
rz(-2.383411) q[0];
sx q[0];
rz(-0.36769205) q[0];
sx q[0];
rz(-1.0962857) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(-0.23263045) q[1];
sx q[1];
rz(-2.0356324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5231424) q[0];
sx q[0];
rz(-1.179427) q[0];
sx q[0];
rz(-0.45446332) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6673379) q[2];
sx q[2];
rz(-1.8372922) q[2];
sx q[2];
rz(1.0702373) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3333862) q[1];
sx q[1];
rz(-1.5667586) q[1];
sx q[1];
rz(1.5334849) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1177103) q[3];
sx q[3];
rz(-2.8571354) q[3];
sx q[3];
rz(0.44265714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37630633) q[2];
sx q[2];
rz(-2.6708965) q[2];
sx q[2];
rz(2.7619696) q[2];
rz(-1.6342573) q[3];
sx q[3];
rz(-1.0999271) q[3];
sx q[3];
rz(-2.1716993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64604243) q[0];
sx q[0];
rz(-2.8046785) q[0];
sx q[0];
rz(-2.2406793) q[0];
rz(-3.0100929) q[1];
sx q[1];
rz(-1.9410746) q[1];
sx q[1];
rz(-0.44201717) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7666047) q[0];
sx q[0];
rz(-0.33549136) q[0];
sx q[0];
rz(-0.25095244) q[0];
x q[1];
rz(2.9346048) q[2];
sx q[2];
rz(-0.0047193165) q[2];
sx q[2];
rz(-2.8568792) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9513248) q[1];
sx q[1];
rz(-0.76808483) q[1];
sx q[1];
rz(-0.70835857) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6824746) q[3];
sx q[3];
rz(-1.1726716) q[3];
sx q[3];
rz(-2.3071837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3758731) q[2];
sx q[2];
rz(-1.7123545) q[2];
sx q[2];
rz(2.5326552) q[2];
rz(-2.7385312) q[3];
sx q[3];
rz(-1.9654704) q[3];
sx q[3];
rz(-2.6398931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.32560638) q[0];
sx q[0];
rz(-2.7222962) q[0];
sx q[0];
rz(2.0035279) q[0];
rz(-1.552938) q[1];
sx q[1];
rz(-2.373003) q[1];
sx q[1];
rz(-0.80702153) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8391221) q[0];
sx q[0];
rz(-0.10286504) q[0];
sx q[0];
rz(-2.8737322) q[0];
x q[1];
rz(1.4078232) q[2];
sx q[2];
rz(-2.0103243) q[2];
sx q[2];
rz(-1.7874315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6389966) q[1];
sx q[1];
rz(-1.4514873) q[1];
sx q[1];
rz(2.9519777) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93385796) q[3];
sx q[3];
rz(-2.6818706) q[3];
sx q[3];
rz(0.62593725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1300065) q[2];
sx q[2];
rz(-0.61379543) q[2];
sx q[2];
rz(-1.2137132) q[2];
rz(3.0774434) q[3];
sx q[3];
rz(-1.6545273) q[3];
sx q[3];
rz(-0.00042644342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5478058) q[0];
sx q[0];
rz(-1.8812027) q[0];
sx q[0];
rz(2.2547145) q[0];
rz(0.15874323) q[1];
sx q[1];
rz(-2.6491149) q[1];
sx q[1];
rz(1.278272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87065164) q[0];
sx q[0];
rz(-2.2566911) q[0];
sx q[0];
rz(-2.565388) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5910225) q[2];
sx q[2];
rz(-1.1614369) q[2];
sx q[2];
rz(2.0602202) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.63250971) q[1];
sx q[1];
rz(-0.69493587) q[1];
sx q[1];
rz(1.664723) q[1];
rz(-pi) q[2];
rz(1.015993) q[3];
sx q[3];
rz(-1.699506) q[3];
sx q[3];
rz(0.57818613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0025803) q[2];
sx q[2];
rz(-2.2860892) q[2];
sx q[2];
rz(2.1878237) q[2];
rz(-1.9468797) q[3];
sx q[3];
rz(-0.84269968) q[3];
sx q[3];
rz(0.4755303) q[3];
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
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82763982) q[0];
sx q[0];
rz(-0.71641818) q[0];
sx q[0];
rz(-2.5415976) q[0];
rz(-2.4586239) q[1];
sx q[1];
rz(-2.3536847) q[1];
sx q[1];
rz(2.8111828) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70008695) q[0];
sx q[0];
rz(-1.7528025) q[0];
sx q[0];
rz(0.84199961) q[0];
x q[1];
rz(-2.1938964) q[2];
sx q[2];
rz(-2.8665339) q[2];
sx q[2];
rz(-0.81240053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.85530419) q[1];
sx q[1];
rz(-1.0897314) q[1];
sx q[1];
rz(0.21800133) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86946671) q[3];
sx q[3];
rz(-1.895213) q[3];
sx q[3];
rz(2.3611435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43634513) q[2];
sx q[2];
rz(-1.0785495) q[2];
sx q[2];
rz(2.516563) q[2];
rz(-2.446512) q[3];
sx q[3];
rz(-1.2806226) q[3];
sx q[3];
rz(3.0888016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2136252) q[0];
sx q[0];
rz(-1.9897505) q[0];
sx q[0];
rz(1.7684162) q[0];
rz(-1.7534509) q[1];
sx q[1];
rz(-0.87166798) q[1];
sx q[1];
rz(-1.6067827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2508188) q[0];
sx q[0];
rz(-1.6060595) q[0];
sx q[0];
rz(1.4232487) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8839726) q[2];
sx q[2];
rz(-1.3968236) q[2];
sx q[2];
rz(-0.18944511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2590547) q[1];
sx q[1];
rz(-1.1657622) q[1];
sx q[1];
rz(-0.71230131) q[1];
rz(-1.2535048) q[3];
sx q[3];
rz(-0.41966885) q[3];
sx q[3];
rz(-0.19007006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9390949) q[2];
sx q[2];
rz(-1.9715344) q[2];
sx q[2];
rz(1.6004174) q[2];
rz(-2.1206858) q[3];
sx q[3];
rz(-1.7424135) q[3];
sx q[3];
rz(-2.8405564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0114667) q[0];
sx q[0];
rz(-0.13700329) q[0];
sx q[0];
rz(-0.89114183) q[0];
rz(-2.4961684) q[1];
sx q[1];
rz(-1.7452469) q[1];
sx q[1];
rz(-2.3088764) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.010983) q[0];
sx q[0];
rz(-2.0314706) q[0];
sx q[0];
rz(1.4185262) q[0];
rz(-pi) q[1];
rz(-1.662976) q[2];
sx q[2];
rz(-2.5155009) q[2];
sx q[2];
rz(-0.41958671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.558185) q[1];
sx q[1];
rz(-1.3761569) q[1];
sx q[1];
rz(2.6411112) q[1];
rz(-pi) q[2];
rz(2.6012558) q[3];
sx q[3];
rz(-1.9129921) q[3];
sx q[3];
rz(-1.0409174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.58925313) q[2];
sx q[2];
rz(-0.66698843) q[2];
sx q[2];
rz(1.5480631) q[2];
rz(0.42406905) q[3];
sx q[3];
rz(-1.4196906) q[3];
sx q[3];
rz(-0.29485318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1770723) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(-0.82343423) q[0];
rz(-1.9442762) q[1];
sx q[1];
rz(-1.4875393) q[1];
sx q[1];
rz(-1.6808602) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4890503) q[0];
sx q[0];
rz(-2.8038137) q[0];
sx q[0];
rz(-1.9073639) q[0];
rz(-1.3235739) q[2];
sx q[2];
rz(-0.6073911) q[2];
sx q[2];
rz(-0.045890778) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9842661) q[1];
sx q[1];
rz(-1.6132024) q[1];
sx q[1];
rz(0.66161109) q[1];
x q[2];
rz(-0.36007215) q[3];
sx q[3];
rz(-1.2692034) q[3];
sx q[3];
rz(-2.2693279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0921649) q[2];
sx q[2];
rz(-0.66114134) q[2];
sx q[2];
rz(3.0180422) q[2];
rz(2.3030247) q[3];
sx q[3];
rz(-2.0288012) q[3];
sx q[3];
rz(-2.6113094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11414828) q[0];
sx q[0];
rz(-1.0076948) q[0];
sx q[0];
rz(1.7281519) q[0];
rz(2.24276) q[1];
sx q[1];
rz(-2.2347968) q[1];
sx q[1];
rz(-2.5487505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44594793) q[0];
sx q[0];
rz(-0.88291016) q[0];
sx q[0];
rz(2.3700729) q[0];
x q[1];
rz(-2.7024621) q[2];
sx q[2];
rz(-0.20072099) q[2];
sx q[2];
rz(1.8619271) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71708116) q[1];
sx q[1];
rz(-0.53679481) q[1];
sx q[1];
rz(2.2473141) q[1];
rz(2.1778706) q[3];
sx q[3];
rz(-2.1603909) q[3];
sx q[3];
rz(2.251925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2841407) q[2];
sx q[2];
rz(-2.393674) q[2];
sx q[2];
rz(0.22135529) q[2];
rz(1.0434693) q[3];
sx q[3];
rz(-2.2011493) q[3];
sx q[3];
rz(-2.7531457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2986044) q[0];
sx q[0];
rz(-2.8589111) q[0];
sx q[0];
rz(-0.84841949) q[0];
rz(0.09659718) q[1];
sx q[1];
rz(-1.074147) q[1];
sx q[1];
rz(-0.75327795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3101075) q[0];
sx q[0];
rz(-1.3436755) q[0];
sx q[0];
rz(-0.102553) q[0];
rz(-3.1084849) q[2];
sx q[2];
rz(-0.82947846) q[2];
sx q[2];
rz(-2.884609) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92434525) q[1];
sx q[1];
rz(-1.0848404) q[1];
sx q[1];
rz(-2.7107581) q[1];
x q[2];
rz(2.4270699) q[3];
sx q[3];
rz(-1.9252637) q[3];
sx q[3];
rz(-2.9936341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0640556) q[2];
sx q[2];
rz(-2.3040743) q[2];
sx q[2];
rz(-0.45292863) q[2];
rz(2.5405267) q[3];
sx q[3];
rz(-1.6744924) q[3];
sx q[3];
rz(-0.99630228) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8375028) q[0];
sx q[0];
rz(-1.4488198) q[0];
sx q[0];
rz(-2.6813843) q[0];
rz(-0.17644633) q[1];
sx q[1];
rz(-2.9063168) q[1];
sx q[1];
rz(1.6520687) q[1];
rz(0.0017133102) q[2];
sx q[2];
rz(-2.7929577) q[2];
sx q[2];
rz(-2.6996573) q[2];
rz(-2.4963958) q[3];
sx q[3];
rz(-2.3020036) q[3];
sx q[3];
rz(0.19097542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

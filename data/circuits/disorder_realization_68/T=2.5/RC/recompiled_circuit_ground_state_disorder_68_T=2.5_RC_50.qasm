OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3026128) q[0];
sx q[0];
rz(-1.4361199) q[0];
sx q[0];
rz(0.21685313) q[0];
rz(2.6537553) q[1];
sx q[1];
rz(-2.2626329) q[1];
sx q[1];
rz(0.15329696) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3756806) q[0];
sx q[0];
rz(-1.3049676) q[0];
sx q[0];
rz(-1.8701118) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6716483) q[2];
sx q[2];
rz(-1.1026088) q[2];
sx q[2];
rz(-0.17344698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5342321) q[1];
sx q[1];
rz(-0.56391729) q[1];
sx q[1];
rz(2.6644071) q[1];
rz(-pi) q[2];
rz(0.80634597) q[3];
sx q[3];
rz(-1.4571179) q[3];
sx q[3];
rz(1.0858842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6903901) q[2];
sx q[2];
rz(-2.5610552) q[2];
sx q[2];
rz(-0.85980493) q[2];
rz(2.9016923) q[3];
sx q[3];
rz(-0.71687117) q[3];
sx q[3];
rz(-0.2828323) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4848223) q[0];
sx q[0];
rz(-1.7253933) q[0];
sx q[0];
rz(-2.2221478) q[0];
rz(2.2942309) q[1];
sx q[1];
rz(-0.58440009) q[1];
sx q[1];
rz(1.1712317) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9860476) q[0];
sx q[0];
rz(-1.520507) q[0];
sx q[0];
rz(1.8296479) q[0];
x q[1];
rz(-1.9483826) q[2];
sx q[2];
rz(-2.3899357) q[2];
sx q[2];
rz(-2.653459) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6714496) q[1];
sx q[1];
rz(-0.43755119) q[1];
sx q[1];
rz(2.2145443) q[1];
rz(-pi) q[2];
x q[2];
rz(0.036072091) q[3];
sx q[3];
rz(-2.6031446) q[3];
sx q[3];
rz(-0.014156646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0688811) q[2];
sx q[2];
rz(-0.88901192) q[2];
sx q[2];
rz(1.2129126) q[2];
rz(-0.21240258) q[3];
sx q[3];
rz(-1.1143149) q[3];
sx q[3];
rz(2.4685278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251959) q[0];
sx q[0];
rz(-1.238751) q[0];
sx q[0];
rz(-2.5585001) q[0];
rz(1.8816226) q[1];
sx q[1];
rz(-2.5446353) q[1];
sx q[1];
rz(-1.3139542) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77480097) q[0];
sx q[0];
rz(-1.2108902) q[0];
sx q[0];
rz(-0.33862305) q[0];
rz(-pi) q[1];
x q[1];
rz(2.407786) q[2];
sx q[2];
rz(-1.5467522) q[2];
sx q[2];
rz(0.95685416) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3594234) q[1];
sx q[1];
rz(-0.12453989) q[1];
sx q[1];
rz(-1.729639) q[1];
x q[2];
rz(-1.9989955) q[3];
sx q[3];
rz(-1.4185662) q[3];
sx q[3];
rz(-0.58109944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49695697) q[2];
sx q[2];
rz(-2.3778215) q[2];
sx q[2];
rz(0.53488564) q[2];
rz(-0.48318091) q[3];
sx q[3];
rz(-1.6202241) q[3];
sx q[3];
rz(3.0218637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508761) q[0];
sx q[0];
rz(-0.058598761) q[0];
sx q[0];
rz(1.6706049) q[0];
rz(1.2141256) q[1];
sx q[1];
rz(-1.7045226) q[1];
sx q[1];
rz(2.6953183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69180079) q[0];
sx q[0];
rz(-1.746382) q[0];
sx q[0];
rz(3.0719766) q[0];
rz(-0.52881119) q[2];
sx q[2];
rz(-1.5870278) q[2];
sx q[2];
rz(2.0179786) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7301) q[1];
sx q[1];
rz(-0.92918438) q[1];
sx q[1];
rz(1.8909251) q[1];
rz(-1.5851444) q[3];
sx q[3];
rz(-1.8135462) q[3];
sx q[3];
rz(1.2402463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2176167) q[2];
sx q[2];
rz(-0.46102229) q[2];
sx q[2];
rz(-0.32611845) q[2];
rz(-2.7255132) q[3];
sx q[3];
rz(-1.6385498) q[3];
sx q[3];
rz(0.7777586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.07051) q[0];
sx q[0];
rz(-1.5776881) q[0];
sx q[0];
rz(0.34410205) q[0];
rz(-2.7857419) q[1];
sx q[1];
rz(-0.75332037) q[1];
sx q[1];
rz(0.25486249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0094285) q[0];
sx q[0];
rz(-2.4309506) q[0];
sx q[0];
rz(-0.406213) q[0];
x q[1];
rz(-1.7698257) q[2];
sx q[2];
rz(-2.6712199) q[2];
sx q[2];
rz(0.34188893) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6333511) q[1];
sx q[1];
rz(-0.76741793) q[1];
sx q[1];
rz(2.2240646) q[1];
rz(0.76538442) q[3];
sx q[3];
rz(-1.2684115) q[3];
sx q[3];
rz(-2.9655448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7039589) q[2];
sx q[2];
rz(-2.2767229) q[2];
sx q[2];
rz(-0.09058365) q[2];
rz(0.80638742) q[3];
sx q[3];
rz(-1.3017637) q[3];
sx q[3];
rz(-1.8346132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.059747132) q[0];
sx q[0];
rz(-1.2223926) q[0];
sx q[0];
rz(-0.61862373) q[0];
rz(2.5289358) q[1];
sx q[1];
rz(-0.63059348) q[1];
sx q[1];
rz(2.9887066) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3090749) q[0];
sx q[0];
rz(-0.33402696) q[0];
sx q[0];
rz(2.3502653) q[0];
x q[1];
rz(2.5618951) q[2];
sx q[2];
rz(-1.6561469) q[2];
sx q[2];
rz(0.63324636) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1807993) q[1];
sx q[1];
rz(-1.5174148) q[1];
sx q[1];
rz(-1.175715) q[1];
x q[2];
rz(-3.116947) q[3];
sx q[3];
rz(-0.99089115) q[3];
sx q[3];
rz(-0.49296311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9360518) q[2];
sx q[2];
rz(-0.80465332) q[2];
sx q[2];
rz(2.0118227) q[2];
rz(-1.0063082) q[3];
sx q[3];
rz(-1.8215424) q[3];
sx q[3];
rz(1.4973076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0928918) q[0];
sx q[0];
rz(-2.0755656) q[0];
sx q[0];
rz(-2.997828) q[0];
rz(0.20214209) q[1];
sx q[1];
rz(-2.6136687) q[1];
sx q[1];
rz(-2.0595097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7162299) q[0];
sx q[0];
rz(-0.83304616) q[0];
sx q[0];
rz(-2.0701755) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80777709) q[2];
sx q[2];
rz(-0.91726724) q[2];
sx q[2];
rz(-1.8612922) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6134167) q[1];
sx q[1];
rz(-2.3966463) q[1];
sx q[1];
rz(-1.2553535) q[1];
rz(0.96636678) q[3];
sx q[3];
rz(-0.92741889) q[3];
sx q[3];
rz(2.9288187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9988592) q[2];
sx q[2];
rz(-1.961668) q[2];
sx q[2];
rz(2.414523) q[2];
rz(1.3149698) q[3];
sx q[3];
rz(-1.246779) q[3];
sx q[3];
rz(-1.6453913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24406782) q[0];
sx q[0];
rz(-1.0777363) q[0];
sx q[0];
rz(-1.9986073) q[0];
rz(-2.9179528) q[1];
sx q[1];
rz(-0.63839212) q[1];
sx q[1];
rz(1.2248096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53608928) q[0];
sx q[0];
rz(-1.8538626) q[0];
sx q[0];
rz(2.8686531) q[0];
x q[1];
rz(1.9838422) q[2];
sx q[2];
rz(-2.2487679) q[2];
sx q[2];
rz(1.7630446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1250028) q[1];
sx q[1];
rz(-1.3666808) q[1];
sx q[1];
rz(0.05929534) q[1];
rz(0.16815987) q[3];
sx q[3];
rz(-0.56795299) q[3];
sx q[3];
rz(1.1105383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58330047) q[2];
sx q[2];
rz(-1.6879098) q[2];
sx q[2];
rz(0.80201403) q[2];
rz(-1.2174886) q[3];
sx q[3];
rz(-3.0018482) q[3];
sx q[3];
rz(1.8961934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32945803) q[0];
sx q[0];
rz(-1.3690925) q[0];
sx q[0];
rz(-1.5611956) q[0];
rz(-0.87248692) q[1];
sx q[1];
rz(-0.28631887) q[1];
sx q[1];
rz(-2.7365541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1096538) q[0];
sx q[0];
rz(-0.73470107) q[0];
sx q[0];
rz(2.3959909) q[0];
rz(1.4680176) q[2];
sx q[2];
rz(-0.98066247) q[2];
sx q[2];
rz(2.0128606) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.26012684) q[1];
sx q[1];
rz(-0.88231371) q[1];
sx q[1];
rz(-2.3860809) q[1];
rz(-1.5902577) q[3];
sx q[3];
rz(-1.7177337) q[3];
sx q[3];
rz(-2.9732239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8672455) q[2];
sx q[2];
rz(-2.2623107) q[2];
sx q[2];
rz(-2.530063) q[2];
rz(1.8379755) q[3];
sx q[3];
rz(-0.36208624) q[3];
sx q[3];
rz(1.1360137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5030454) q[0];
sx q[0];
rz(-1.2940116) q[0];
sx q[0];
rz(3.0882623) q[0];
rz(0.57180014) q[1];
sx q[1];
rz(-1.5170393) q[1];
sx q[1];
rz(1.2791876) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93829556) q[0];
sx q[0];
rz(-0.38559993) q[0];
sx q[0];
rz(2.4176265) q[0];
x q[1];
rz(-1.0903075) q[2];
sx q[2];
rz(-1.4391403) q[2];
sx q[2];
rz(2.5327794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0030583) q[1];
sx q[1];
rz(-1.492625) q[1];
sx q[1];
rz(2.1581236) q[1];
rz(-pi) q[2];
rz(1.199715) q[3];
sx q[3];
rz(-1.7744007) q[3];
sx q[3];
rz(2.9534087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46644396) q[2];
sx q[2];
rz(-1.5326591) q[2];
sx q[2];
rz(-0.68501985) q[2];
rz(-2.3576665) q[3];
sx q[3];
rz(-0.42015606) q[3];
sx q[3];
rz(-2.4043731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53981275) q[0];
sx q[0];
rz(-1.2417326) q[0];
sx q[0];
rz(-1.7057521) q[0];
rz(2.336179) q[1];
sx q[1];
rz(-0.46331159) q[1];
sx q[1];
rz(0.70809271) q[1];
rz(-0.63486256) q[2];
sx q[2];
rz(-0.35463968) q[2];
sx q[2];
rz(-1.126033) q[2];
rz(0.41498923) q[3];
sx q[3];
rz(-1.4616953) q[3];
sx q[3];
rz(-2.2898522) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

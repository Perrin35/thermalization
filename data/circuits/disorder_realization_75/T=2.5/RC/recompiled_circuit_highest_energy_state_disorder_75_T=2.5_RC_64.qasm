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
rz(1.6504352) q[0];
sx q[0];
rz(2.7661134) q[0];
sx q[0];
rz(9.8162415) q[0];
rz(0.23671167) q[1];
sx q[1];
rz(-2.1032636) q[1];
sx q[1];
rz(-2.2957323) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0856871) q[0];
sx q[0];
rz(-2.2330267) q[0];
sx q[0];
rz(-1.7895749) q[0];
x q[1];
rz(1.1089788) q[2];
sx q[2];
rz(-1.9792288) q[2];
sx q[2];
rz(-1.7866745) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7995092) q[1];
sx q[1];
rz(-1.2984283) q[1];
sx q[1];
rz(1.8567159) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99669807) q[3];
sx q[3];
rz(-0.59060589) q[3];
sx q[3];
rz(-2.8656538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5819431) q[2];
sx q[2];
rz(-2.1008284) q[2];
sx q[2];
rz(-0.2156269) q[2];
rz(0.98254472) q[3];
sx q[3];
rz(-1.6590786) q[3];
sx q[3];
rz(1.2561579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243645) q[0];
sx q[0];
rz(-3.0475782) q[0];
sx q[0];
rz(-2.7650058) q[0];
rz(-1.4812034) q[1];
sx q[1];
rz(-1.9638289) q[1];
sx q[1];
rz(-1.0053763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97949147) q[0];
sx q[0];
rz(-2.8164356) q[0];
sx q[0];
rz(-1.5942911) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2475904) q[2];
sx q[2];
rz(-1.4419705) q[2];
sx q[2];
rz(1.027838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0223234) q[1];
sx q[1];
rz(-0.6695153) q[1];
sx q[1];
rz(-3.0864848) q[1];
x q[2];
rz(2.5915089) q[3];
sx q[3];
rz(-1.2445916) q[3];
sx q[3];
rz(-2.4334513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5873831) q[2];
sx q[2];
rz(-2.6770834) q[2];
sx q[2];
rz(-2.0655538) q[2];
rz(2.6739323) q[3];
sx q[3];
rz(-2.3105919) q[3];
sx q[3];
rz(0.20461288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6429546) q[0];
sx q[0];
rz(-1.0774281) q[0];
sx q[0];
rz(-2.0306008) q[0];
rz(-0.30771646) q[1];
sx q[1];
rz(-1.6650763) q[1];
sx q[1];
rz(-1.3002546) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1210403) q[0];
sx q[0];
rz(-2.5835134) q[0];
sx q[0];
rz(0.16772018) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57794364) q[2];
sx q[2];
rz(-1.6010849) q[2];
sx q[2];
rz(-0.89627111) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3798602) q[1];
sx q[1];
rz(-0.76660313) q[1];
sx q[1];
rz(2.2025781) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6407317) q[3];
sx q[3];
rz(-1.259049) q[3];
sx q[3];
rz(2.4846128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51851455) q[2];
sx q[2];
rz(-2.0936421) q[2];
sx q[2];
rz(-2.9497362) q[2];
rz(1.1040556) q[3];
sx q[3];
rz(-2.7147229) q[3];
sx q[3];
rz(1.9328611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0693102) q[0];
sx q[0];
rz(-0.59232124) q[0];
sx q[0];
rz(2.2533672) q[0];
rz(2.8690673) q[1];
sx q[1];
rz(-1.2963632) q[1];
sx q[1];
rz(-1.9857508) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10242045) q[0];
sx q[0];
rz(-1.8646324) q[0];
sx q[0];
rz(-3.1320024) q[0];
rz(-pi) q[1];
rz(-1.5290909) q[2];
sx q[2];
rz(-1.0805849) q[2];
sx q[2];
rz(0.34402564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0371746) q[1];
sx q[1];
rz(-1.2942477) q[1];
sx q[1];
rz(2.3642703) q[1];
rz(1.5748197) q[3];
sx q[3];
rz(-2.0186264) q[3];
sx q[3];
rz(-2.8886917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3027975) q[2];
sx q[2];
rz(-0.28442997) q[2];
sx q[2];
rz(2.3885041) q[2];
rz(2.6436842) q[3];
sx q[3];
rz(-1.6585191) q[3];
sx q[3];
rz(-2.8411617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4860151) q[0];
sx q[0];
rz(-2.1652554) q[0];
sx q[0];
rz(1.2855541) q[0];
rz(1.1898419) q[1];
sx q[1];
rz(-2.4557476) q[1];
sx q[1];
rz(1.5815585) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6356892) q[0];
sx q[0];
rz(-2.990953) q[0];
sx q[0];
rz(-2.1156359) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6597775) q[2];
sx q[2];
rz(-2.8236532) q[2];
sx q[2];
rz(-0.16127333) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1413708) q[1];
sx q[1];
rz(-2.0061778) q[1];
sx q[1];
rz(1.8069827) q[1];
rz(-2.8772379) q[3];
sx q[3];
rz(-1.3761545) q[3];
sx q[3];
rz(2.1771012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2994284) q[2];
sx q[2];
rz(-2.2101768) q[2];
sx q[2];
rz(3.0038339) q[2];
rz(-2.6970741) q[3];
sx q[3];
rz(-1.7701745) q[3];
sx q[3];
rz(-0.81454149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21823068) q[0];
sx q[0];
rz(-2.2683999) q[0];
sx q[0];
rz(-0.55106226) q[0];
rz(-0.29523826) q[1];
sx q[1];
rz(-2.4199838) q[1];
sx q[1];
rz(-0.74337426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0966885) q[0];
sx q[0];
rz(-1.7040692) q[0];
sx q[0];
rz(-3.0267984) q[0];
rz(-0.34140519) q[2];
sx q[2];
rz(-1.3479173) q[2];
sx q[2];
rz(-0.49431942) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8507) q[1];
sx q[1];
rz(-2.2883483) q[1];
sx q[1];
rz(0.13137551) q[1];
rz(-2.4427209) q[3];
sx q[3];
rz(-2.0783608) q[3];
sx q[3];
rz(2.5273049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4117406) q[2];
sx q[2];
rz(-2.9509632) q[2];
sx q[2];
rz(-2.8594678) q[2];
rz(-2.7052346) q[3];
sx q[3];
rz(-1.3151508) q[3];
sx q[3];
rz(-2.8859629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2082763) q[0];
sx q[0];
rz(-2.0774807) q[0];
sx q[0];
rz(-2.5489885) q[0];
rz(1.1366049) q[1];
sx q[1];
rz(-1.2717783) q[1];
sx q[1];
rz(-2.073854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.655433) q[0];
sx q[0];
rz(-2.2256333) q[0];
sx q[0];
rz(1.372712) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4456533) q[2];
sx q[2];
rz(-1.9370902) q[2];
sx q[2];
rz(-0.83331457) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6954036) q[1];
sx q[1];
rz(-2.0430889) q[1];
sx q[1];
rz(-0.31361927) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41530825) q[3];
sx q[3];
rz(-1.7303932) q[3];
sx q[3];
rz(0.94091614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9390255) q[2];
sx q[2];
rz(-1.86684) q[2];
sx q[2];
rz(-2.1620046) q[2];
rz(-3.0136287) q[3];
sx q[3];
rz(-2.1991859) q[3];
sx q[3];
rz(-1.6894107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1487883) q[0];
sx q[0];
rz(-2.1908741) q[0];
sx q[0];
rz(0.08091452) q[0];
rz(-0.12123904) q[1];
sx q[1];
rz(-2.464005) q[1];
sx q[1];
rz(1.1239207) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8351562) q[0];
sx q[0];
rz(-0.84266169) q[0];
sx q[0];
rz(2.6373498) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6155682) q[2];
sx q[2];
rz(-1.977619) q[2];
sx q[2];
rz(-1.3699832) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2595265) q[1];
sx q[1];
rz(-1.9212544) q[1];
sx q[1];
rz(-1.6762074) q[1];
rz(-pi) q[2];
rz(2.9714576) q[3];
sx q[3];
rz(-2.0960582) q[3];
sx q[3];
rz(-0.13290031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27440444) q[2];
sx q[2];
rz(-0.97325456) q[2];
sx q[2];
rz(-1.493206) q[2];
rz(2.6793206) q[3];
sx q[3];
rz(-2.3267764) q[3];
sx q[3];
rz(1.7260684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82861519) q[0];
sx q[0];
rz(-1.7531489) q[0];
sx q[0];
rz(1.9762565) q[0];
rz(3.0997711) q[1];
sx q[1];
rz(-0.98894293) q[1];
sx q[1];
rz(1.8080541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81105574) q[0];
sx q[0];
rz(-1.0687318) q[0];
sx q[0];
rz(-0.83107019) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8227764) q[2];
sx q[2];
rz(-1.5747264) q[2];
sx q[2];
rz(1.9424644) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11631498) q[1];
sx q[1];
rz(-1.9851369) q[1];
sx q[1];
rz(-1.6700909) q[1];
rz(-pi) q[2];
rz(-2.4163212) q[3];
sx q[3];
rz(-2.9915611) q[3];
sx q[3];
rz(1.6440441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8970783) q[2];
sx q[2];
rz(-1.3099193) q[2];
sx q[2];
rz(0.46373996) q[2];
rz(0.63716277) q[3];
sx q[3];
rz(-1.5183247) q[3];
sx q[3];
rz(0.37850982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8265182) q[0];
sx q[0];
rz(-1.0239064) q[0];
sx q[0];
rz(1.6465323) q[0];
rz(0.82978326) q[1];
sx q[1];
rz(-1.8779495) q[1];
sx q[1];
rz(-1.8216546) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0129725) q[0];
sx q[0];
rz(-1.2018328) q[0];
sx q[0];
rz(1.4340543) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7097335) q[2];
sx q[2];
rz(-1.3629902) q[2];
sx q[2];
rz(-0.49655216) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1935008) q[1];
sx q[1];
rz(-1.337916) q[1];
sx q[1];
rz(-1.7786673) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9687555) q[3];
sx q[3];
rz(-2.5423942) q[3];
sx q[3];
rz(-2.9323682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39184555) q[2];
sx q[2];
rz(-1.7228935) q[2];
sx q[2];
rz(-0.5113655) q[2];
rz(0.5558719) q[3];
sx q[3];
rz(-2.0248196) q[3];
sx q[3];
rz(-2.6218124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2340672) q[0];
sx q[0];
rz(-3.0341442) q[0];
sx q[0];
rz(-2.7035614) q[0];
rz(3.0829433) q[1];
sx q[1];
rz(-1.6390683) q[1];
sx q[1];
rz(-1.0674089) q[1];
rz(-1.5127586) q[2];
sx q[2];
rz(-0.87067247) q[2];
sx q[2];
rz(-1.6280328) q[2];
rz(0.41154804) q[3];
sx q[3];
rz(-0.66933142) q[3];
sx q[3];
rz(0.29585024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

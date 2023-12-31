OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(-0.056161031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4958772) q[0];
sx q[0];
rz(-2.0464532) q[0];
sx q[0];
rz(-0.23947421) q[0];
x q[1];
rz(1.2916318) q[2];
sx q[2];
rz(-2.4764875) q[2];
sx q[2];
rz(-1.6601738) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.010477) q[1];
sx q[1];
rz(-1.8382204) q[1];
sx q[1];
rz(-2.1285776) q[1];
rz(-pi) q[2];
rz(1.6151186) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(2.0093902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-0.4326694) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(-2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(1.312785) q[0];
rz(2.9361172) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.9899433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1412927) q[0];
sx q[0];
rz(-0.89501689) q[0];
sx q[0];
rz(-0.38210259) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58354124) q[2];
sx q[2];
rz(-1.984664) q[2];
sx q[2];
rz(-1.5838503) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4558251) q[1];
sx q[1];
rz(-2.0058504) q[1];
sx q[1];
rz(-2.0352092) q[1];
rz(-pi) q[2];
rz(2.6395256) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(-1.7934007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.0097222086) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(3.1047399) q[0];
rz(-0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(-0.056578606) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73233561) q[0];
sx q[0];
rz(-1.1261254) q[0];
sx q[0];
rz(2.1952941) q[0];
rz(-1.4372196) q[2];
sx q[2];
rz(-1.8250416) q[2];
sx q[2];
rz(-2.996252) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0128855) q[1];
sx q[1];
rz(-1.7696847) q[1];
sx q[1];
rz(-0.30602869) q[1];
x q[2];
rz(-1.1832841) q[3];
sx q[3];
rz(-0.96252493) q[3];
sx q[3];
rz(-1.1063948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-2.2154714) q[2];
rz(-2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2264003) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(-1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-2.6779968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47953654) q[0];
sx q[0];
rz(-0.89131309) q[0];
sx q[0];
rz(-0.38328538) q[0];
x q[1];
rz(1.6323339) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(-2.5663944) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0193034) q[1];
sx q[1];
rz(-0.79320723) q[1];
sx q[1];
rz(1.8122458) q[1];
x q[2];
rz(-1.0158402) q[3];
sx q[3];
rz(-1.1557126) q[3];
sx q[3];
rz(-2.2500028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(0.049499361) q[2];
rz(-0.1285304) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(-0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(0.4822576) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(0.94435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34236318) q[0];
sx q[0];
rz(-0.21393299) q[0];
sx q[0];
rz(1.7844723) q[0];
rz(-pi) q[1];
rz(1.1733426) q[2];
sx q[2];
rz(-0.83565088) q[2];
sx q[2];
rz(-3.1189001) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9152865) q[1];
sx q[1];
rz(-1.5686791) q[1];
sx q[1];
rz(1.6318984) q[1];
x q[2];
rz(0.51575757) q[3];
sx q[3];
rz(-2.6532647) q[3];
sx q[3];
rz(2.8749136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.258761) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(3.0333701) q[2];
rz(-3.1392858) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(-0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43679431) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(0.58445245) q[0];
rz(0.8862409) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(0.054919682) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45860896) q[0];
sx q[0];
rz(-2.8603472) q[0];
sx q[0];
rz(-1.8722697) q[0];
rz(-pi) q[1];
rz(0.018718406) q[2];
sx q[2];
rz(-1.8939549) q[2];
sx q[2];
rz(1.2013555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1850486) q[1];
sx q[1];
rz(-1.946432) q[1];
sx q[1];
rz(2.5513785) q[1];
x q[2];
rz(0.46338007) q[3];
sx q[3];
rz(-1.0372835) q[3];
sx q[3];
rz(2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(2.1248655) q[2];
rz(0.54404849) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(2.8570535) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-2.231266) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0383354) q[0];
sx q[0];
rz(-3.0568125) q[0];
sx q[0];
rz(-0.8707365) q[0];
rz(-pi) q[1];
rz(0.89332135) q[2];
sx q[2];
rz(-2.6258694) q[2];
sx q[2];
rz(0.90781462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42156223) q[1];
sx q[1];
rz(-1.7517462) q[1];
sx q[1];
rz(2.6471495) q[1];
x q[2];
rz(2.4207553) q[3];
sx q[3];
rz(-2.0257054) q[3];
sx q[3];
rz(-2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(0.92203036) q[2];
rz(-2.5743124) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(-0.63240504) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(2.8410889) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757272) q[0];
sx q[0];
rz(-2.2458796) q[0];
sx q[0];
rz(0.48203326) q[0];
x q[1];
rz(0.6408765) q[2];
sx q[2];
rz(-2.2705728) q[2];
sx q[2];
rz(-1.7388294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.955519) q[1];
sx q[1];
rz(-1.7454073) q[1];
sx q[1];
rz(-0.34592918) q[1];
rz(-pi) q[2];
x q[2];
rz(2.825533) q[3];
sx q[3];
rz(-0.83871597) q[3];
sx q[3];
rz(0.72533208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.58632103) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(0.78197455) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(0.12776275) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-0.75884563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42620537) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(-0.018652648) q[0];
rz(-pi) q[1];
rz(-1.7636289) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(2.6892975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89275937) q[1];
sx q[1];
rz(-1.7264688) q[1];
sx q[1];
rz(2.4397736) q[1];
rz(-pi) q[2];
rz(-0.90518732) q[3];
sx q[3];
rz(-1.5822516) q[3];
sx q[3];
rz(-0.36597914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(-0.49003595) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(3.066257) q[0];
rz(0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-0.60992253) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41994914) q[0];
sx q[0];
rz(-2.3099265) q[0];
sx q[0];
rz(-1.182611) q[0];
x q[1];
rz(0.37656017) q[2];
sx q[2];
rz(-1.3137523) q[2];
sx q[2];
rz(-0.10274796) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8049106) q[1];
sx q[1];
rz(-1.33178) q[1];
sx q[1];
rz(-2.3906624) q[1];
rz(0.53819733) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.909409) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(-0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(0.65080416) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.4177633) q[2];
sx q[2];
rz(-0.19822181) q[2];
sx q[2];
rz(-0.76186686) q[2];
rz(0.25272947) q[3];
sx q[3];
rz(-1.1128294) q[3];
sx q[3];
rz(2.2313234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

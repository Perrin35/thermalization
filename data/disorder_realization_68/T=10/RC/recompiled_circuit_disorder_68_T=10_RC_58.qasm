OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(-0.49180254) q[0];
sx q[0];
rz(0.1879745) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(-1.6245276) q[1];
sx q[1];
rz(-2.7741073) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7029019) q[0];
sx q[0];
rz(-2.5950948) q[0];
sx q[0];
rz(-1.7706857) q[0];
x q[1];
rz(-1.3221402) q[2];
sx q[2];
rz(-2.6373632) q[2];
sx q[2];
rz(2.3766975) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.623466) q[1];
sx q[1];
rz(-1.4033485) q[1];
sx q[1];
rz(1.3256339) q[1];
rz(1.3931307) q[3];
sx q[3];
rz(-1.2347617) q[3];
sx q[3];
rz(-2.6538268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.964103) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(1.8356813) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630163) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(-2.6696894) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(-2.205251) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72915709) q[0];
sx q[0];
rz(-1.658174) q[0];
sx q[0];
rz(-0.23075128) q[0];
rz(-2.721644) q[2];
sx q[2];
rz(-2.311085) q[2];
sx q[2];
rz(0.74479693) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.7649819) q[1];
sx q[1];
rz(-2.2602343) q[1];
sx q[1];
rz(1.7060243) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0735047) q[3];
sx q[3];
rz(-1.5912676) q[3];
sx q[3];
rz(-0.81077829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77461809) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(-1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(-0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24580978) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(-2.202503) q[0];
rz(-0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(2.5476707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.827841) q[0];
sx q[0];
rz(-1.272164) q[0];
sx q[0];
rz(-1.2350425) q[0];
rz(2.2592696) q[2];
sx q[2];
rz(-1.5261298) q[2];
sx q[2];
rz(2.4633212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8507104) q[1];
sx q[1];
rz(-1.9899273) q[1];
sx q[1];
rz(-0.23371975) q[1];
x q[2];
rz(2.6731554) q[3];
sx q[3];
rz(-1.2611946) q[3];
sx q[3];
rz(2.7999511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64017355) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(1.7017986) q[2];
rz(-2.7539608) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(-2.7534527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.7664292) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(2.638812) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(-0.75685135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3319791) q[0];
sx q[0];
rz(-0.50052128) q[0];
sx q[0];
rz(-2.39397) q[0];
rz(-pi) q[1];
rz(-1.5449764) q[2];
sx q[2];
rz(-0.46586793) q[2];
sx q[2];
rz(-0.70170882) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5731569) q[1];
sx q[1];
rz(-2.4385298) q[1];
sx q[1];
rz(2.1431124) q[1];
rz(0.54550708) q[3];
sx q[3];
rz(-2.3295998) q[3];
sx q[3];
rz(3.0363887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42671529) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(1.4871917) q[2];
rz(-2.5590844) q[3];
sx q[3];
rz(-1.094386) q[3];
sx q[3];
rz(0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-2.0563545) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(1.2879397) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(-2.0910738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8331063) q[0];
sx q[0];
rz(-1.3163438) q[0];
sx q[0];
rz(1.2480877) q[0];
rz(-pi) q[1];
rz(0.44287037) q[2];
sx q[2];
rz(-2.1211229) q[2];
sx q[2];
rz(1.0843104) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1615636) q[1];
sx q[1];
rz(-2.1254351) q[1];
sx q[1];
rz(0.6273004) q[1];
x q[2];
rz(0.80661185) q[3];
sx q[3];
rz(-1.9730554) q[3];
sx q[3];
rz(0.81392399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.918255) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(0.63344947) q[2];
rz(-1.194362) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.69960064) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(-1.6794499) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1519449) q[0];
sx q[0];
rz(-1.6501353) q[0];
sx q[0];
rz(-0.21110714) q[0];
rz(-pi) q[1];
rz(2.2199549) q[2];
sx q[2];
rz(-1.7643133) q[2];
sx q[2];
rz(2.9606539) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.51965442) q[1];
sx q[1];
rz(-0.42236537) q[1];
sx q[1];
rz(0.42610355) q[1];
x q[2];
rz(-1.2268279) q[3];
sx q[3];
rz(-1.868639) q[3];
sx q[3];
rz(2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9399461) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(2.3366826) q[2];
rz(-1.9645875) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(1.0837519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8686304) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(-0.92765635) q[0];
rz(-2.1169128) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(-1.0120846) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18400684) q[0];
sx q[0];
rz(-2.036096) q[0];
sx q[0];
rz(2.1561949) q[0];
rz(-pi) q[1];
rz(2.8132642) q[2];
sx q[2];
rz(-2.7331181) q[2];
sx q[2];
rz(-2.7834053) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95841366) q[1];
sx q[1];
rz(-1.603754) q[1];
sx q[1];
rz(-2.5977913) q[1];
x q[2];
rz(-1.6530232) q[3];
sx q[3];
rz(-1.0136908) q[3];
sx q[3];
rz(1.132387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6107789) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(-2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(-2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5291418) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(0.21417831) q[0];
rz(-2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(-2.8578551) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8026233) q[0];
sx q[0];
rz(-1.5366652) q[0];
sx q[0];
rz(-0.3011093) q[0];
rz(-pi) q[1];
rz(2.0869414) q[2];
sx q[2];
rz(-0.37545855) q[2];
sx q[2];
rz(1.6106538) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1188018) q[1];
sx q[1];
rz(-2.6036501) q[1];
sx q[1];
rz(-2.7522037) q[1];
rz(-0.92054263) q[3];
sx q[3];
rz(-1.0245819) q[3];
sx q[3];
rz(-1.2342681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6909137) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(-1.696375) q[2];
rz(-1.5971659) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.504869) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(0.069256393) q[0];
rz(1.6537369) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(-1.5725296) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3587787) q[0];
sx q[0];
rz(-0.94952119) q[0];
sx q[0];
rz(0.57399477) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8157418) q[2];
sx q[2];
rz(-1.4928276) q[2];
sx q[2];
rz(-2.1577912) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.442765) q[1];
sx q[1];
rz(-0.24833939) q[1];
sx q[1];
rz(0.34766867) q[1];
rz(-0.92215718) q[3];
sx q[3];
rz(-1.7382009) q[3];
sx q[3];
rz(2.5006014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9562324) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(3.0017079) q[2];
rz(0.36758962) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-0.39127025) q[0];
sx q[0];
rz(2.4998253) q[0];
rz(-1.9104674) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(2.8737601) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56775996) q[0];
sx q[0];
rz(-1.1288252) q[0];
sx q[0];
rz(2.3964336) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7540625) q[2];
sx q[2];
rz(-1.4472618) q[2];
sx q[2];
rz(1.087041) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.47083449) q[1];
sx q[1];
rz(-1.2467614) q[1];
sx q[1];
rz(1.8761937) q[1];
rz(-pi) q[2];
rz(-2.0249428) q[3];
sx q[3];
rz(-2.5513259) q[3];
sx q[3];
rz(-0.86436194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(-1.4769185) q[2];
rz(-2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.01263604) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(-2.3616882) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(2.4308464) q[2];
sx q[2];
rz(-2.032861) q[2];
sx q[2];
rz(-2.1326333) q[2];
rz(0.6775425) q[3];
sx q[3];
rz(-1.4031706) q[3];
sx q[3];
rz(-1.8085898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

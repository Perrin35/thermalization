OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.087091669) q[0];
sx q[0];
rz(-2.2890685) q[0];
sx q[0];
rz(2.8791715) q[0];
rz(-3.4497058) q[1];
sx q[1];
rz(0.83642712) q[1];
sx q[1];
rz(15.717419) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.165929) q[0];
sx q[0];
rz(-1.0895415) q[0];
sx q[0];
rz(-0.4939618) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0980564) q[2];
sx q[2];
rz(-1.6400717) q[2];
sx q[2];
rz(-0.89292282) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.63543273) q[1];
sx q[1];
rz(-0.80997688) q[1];
sx q[1];
rz(-0.86829888) q[1];
rz(1.057714) q[3];
sx q[3];
rz(-0.82699593) q[3];
sx q[3];
rz(-1.8216474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0414163) q[2];
sx q[2];
rz(-2.2569423) q[2];
sx q[2];
rz(0.46119383) q[2];
rz(2.6395116) q[3];
sx q[3];
rz(-0.82379782) q[3];
sx q[3];
rz(-1.7320777) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13040386) q[0];
sx q[0];
rz(-1.4606425) q[0];
sx q[0];
rz(-2.9050264) q[0];
rz(0.81831167) q[1];
sx q[1];
rz(-1.9383483) q[1];
sx q[1];
rz(-0.94258211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494076) q[0];
sx q[0];
rz(-1.5680917) q[0];
sx q[0];
rz(1.6016141) q[0];
x q[1];
rz(-0.06748345) q[2];
sx q[2];
rz(-1.4500256) q[2];
sx q[2];
rz(-1.6492594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2494804) q[1];
sx q[1];
rz(-1.5722407) q[1];
sx q[1];
rz(2.4066928) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0560741) q[3];
sx q[3];
rz(-1.8746146) q[3];
sx q[3];
rz(-1.3930381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4054823) q[2];
sx q[2];
rz(-2.5030899) q[2];
sx q[2];
rz(-0.68518266) q[2];
rz(-2.3403366) q[3];
sx q[3];
rz(-1.6359436) q[3];
sx q[3];
rz(1.8906458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41246978) q[0];
sx q[0];
rz(-2.6710489) q[0];
sx q[0];
rz(0.98010081) q[0];
rz(0.12067548) q[1];
sx q[1];
rz(-1.9735034) q[1];
sx q[1];
rz(1.6623704) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50608327) q[0];
sx q[0];
rz(-1.4091307) q[0];
sx q[0];
rz(0.76753214) q[0];
rz(2.6363716) q[2];
sx q[2];
rz(-1.1596934) q[2];
sx q[2];
rz(-1.1875325) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2328636) q[1];
sx q[1];
rz(-1.5655978) q[1];
sx q[1];
rz(-0.0089110891) q[1];
rz(-pi) q[2];
rz(-0.5898409) q[3];
sx q[3];
rz(-1.2882659) q[3];
sx q[3];
rz(-0.88385201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7145308) q[2];
sx q[2];
rz(-1.6468628) q[2];
sx q[2];
rz(-2.96116) q[2];
rz(2.7212972) q[3];
sx q[3];
rz(-2.1642978) q[3];
sx q[3];
rz(-0.54496566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412398) q[0];
sx q[0];
rz(-1.350133) q[0];
sx q[0];
rz(0.049276503) q[0];
rz(-0.73206466) q[1];
sx q[1];
rz(-2.2923636) q[1];
sx q[1];
rz(-1.3759618) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7017562) q[0];
sx q[0];
rz(-1.6876093) q[0];
sx q[0];
rz(-2.9573836) q[0];
rz(-1.3591197) q[2];
sx q[2];
rz(-0.81640252) q[2];
sx q[2];
rz(-0.27118452) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.84280864) q[1];
sx q[1];
rz(-2.7453303) q[1];
sx q[1];
rz(2.9664459) q[1];
rz(0.84305544) q[3];
sx q[3];
rz(-0.84765654) q[3];
sx q[3];
rz(1.8384733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10776082) q[2];
sx q[2];
rz(-1.7511448) q[2];
sx q[2];
rz(2.4793009) q[2];
rz(-1.8108588) q[3];
sx q[3];
rz(-2.3190506) q[3];
sx q[3];
rz(1.8432157) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0996284) q[0];
sx q[0];
rz(-0.30655107) q[0];
sx q[0];
rz(2.1339259) q[0];
rz(3.0362466) q[1];
sx q[1];
rz(-2.4844929) q[1];
sx q[1];
rz(-0.23060051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67543519) q[0];
sx q[0];
rz(-0.49336067) q[0];
sx q[0];
rz(2.3301279) q[0];
x q[1];
rz(-0.23955524) q[2];
sx q[2];
rz(-1.2192246) q[2];
sx q[2];
rz(0.4229751) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8457875) q[1];
sx q[1];
rz(-1.6562067) q[1];
sx q[1];
rz(-0.16525646) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19959656) q[3];
sx q[3];
rz(-0.65465876) q[3];
sx q[3];
rz(2.9195291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0577804) q[2];
sx q[2];
rz(-2.1746706) q[2];
sx q[2];
rz(2.9528565) q[2];
rz(1.905929) q[3];
sx q[3];
rz(-0.50714791) q[3];
sx q[3];
rz(1.2602826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96596232) q[0];
sx q[0];
rz(-1.9251134) q[0];
sx q[0];
rz(2.8438582) q[0];
rz(3.0689902) q[1];
sx q[1];
rz(-2.4563792) q[1];
sx q[1];
rz(-2.4593478) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84083623) q[0];
sx q[0];
rz(-1.8634691) q[0];
sx q[0];
rz(-0.93772823) q[0];
rz(1.7715681) q[2];
sx q[2];
rz(-2.4983495) q[2];
sx q[2];
rz(-2.7386576) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.025921019) q[1];
sx q[1];
rz(-2.762315) q[1];
sx q[1];
rz(1.9801407) q[1];
rz(-pi) q[2];
rz(2.7704084) q[3];
sx q[3];
rz(-1.0044349) q[3];
sx q[3];
rz(2.0831828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5630774) q[2];
sx q[2];
rz(-0.35334057) q[2];
sx q[2];
rz(-0.087470857) q[2];
rz(-2.2443917) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(2.7520666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(3.0422269) q[0];
sx q[0];
rz(-1.1023738) q[0];
sx q[0];
rz(-1.1627831) q[0];
rz(2.9258264) q[1];
sx q[1];
rz(-2.3240604) q[1];
sx q[1];
rz(-2.20631) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1638086) q[0];
sx q[0];
rz(-0.80923128) q[0];
sx q[0];
rz(-2.8092395) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0832418) q[2];
sx q[2];
rz(-0.29428665) q[2];
sx q[2];
rz(-2.2541915) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7120119) q[1];
sx q[1];
rz(-1.64974) q[1];
sx q[1];
rz(2.354875) q[1];
rz(0.31511001) q[3];
sx q[3];
rz(-1.2602046) q[3];
sx q[3];
rz(2.9950708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0995348) q[2];
sx q[2];
rz(-2.4877986) q[2];
sx q[2];
rz(-0.81134861) q[2];
rz(2.8387496) q[3];
sx q[3];
rz(-1.9767913) q[3];
sx q[3];
rz(-2.5109049) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52043668) q[0];
sx q[0];
rz(-2.8148837) q[0];
sx q[0];
rz(-2.4503571) q[0];
rz(-0.65903819) q[1];
sx q[1];
rz(-1.6233416) q[1];
sx q[1];
rz(-0.59115994) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1189445) q[0];
sx q[0];
rz(-1.582358) q[0];
sx q[0];
rz(-1.287879) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0760087) q[2];
sx q[2];
rz(-0.55375868) q[2];
sx q[2];
rz(-0.29386917) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5308705) q[1];
sx q[1];
rz(-1.8962144) q[1];
sx q[1];
rz(-2.3476092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.05686111) q[3];
sx q[3];
rz(-1.8627852) q[3];
sx q[3];
rz(2.9546471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14472321) q[2];
sx q[2];
rz(-0.99088061) q[2];
sx q[2];
rz(1.4608176) q[2];
rz(0.75374323) q[3];
sx q[3];
rz(-1.4494579) q[3];
sx q[3];
rz(2.6678705) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4562456) q[0];
sx q[0];
rz(-1.101475) q[0];
sx q[0];
rz(2.9154678) q[0];
rz(-1.6926758) q[1];
sx q[1];
rz(-2.1782404) q[1];
sx q[1];
rz(-2.1400145) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8346092) q[0];
sx q[0];
rz(-0.98794565) q[0];
sx q[0];
rz(-1.7442763) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44160812) q[2];
sx q[2];
rz(-0.749365) q[2];
sx q[2];
rz(-2.6059849) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68943095) q[1];
sx q[1];
rz(-0.2161444) q[1];
sx q[1];
rz(-0.5139022) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4572152) q[3];
sx q[3];
rz(-1.805814) q[3];
sx q[3];
rz(-0.37684611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5292624) q[2];
sx q[2];
rz(-0.77028209) q[2];
sx q[2];
rz(0.15753499) q[2];
rz(-0.15466386) q[3];
sx q[3];
rz(-1.1636584) q[3];
sx q[3];
rz(-2.7111588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73151076) q[0];
sx q[0];
rz(-1.9022576) q[0];
sx q[0];
rz(-2.4391158) q[0];
rz(0.53600535) q[1];
sx q[1];
rz(-2.3119226) q[1];
sx q[1];
rz(1.3943256) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047612543) q[0];
sx q[0];
rz(-1.5965921) q[0];
sx q[0];
rz(1.4071147) q[0];
rz(-pi) q[1];
rz(3.0476116) q[2];
sx q[2];
rz(-1.4620442) q[2];
sx q[2];
rz(0.55545872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6977384) q[1];
sx q[1];
rz(-1.8540253) q[1];
sx q[1];
rz(-2.1069788) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1411885) q[3];
sx q[3];
rz(-2.7260216) q[3];
sx q[3];
rz(-0.11041752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.006762) q[2];
sx q[2];
rz(-2.4595021) q[2];
sx q[2];
rz(-1.6949867) q[2];
rz(2.8805978) q[3];
sx q[3];
rz(-1.010453) q[3];
sx q[3];
rz(-1.235777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21988729) q[0];
sx q[0];
rz(-2.4867262) q[0];
sx q[0];
rz(2.126271) q[0];
rz(-2.0657397) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(0.27412065) q[2];
sx q[2];
rz(-0.39015301) q[2];
sx q[2];
rz(2.9377666) q[2];
rz(0.45810926) q[3];
sx q[3];
rz(-0.81732133) q[3];
sx q[3];
rz(2.2504239) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

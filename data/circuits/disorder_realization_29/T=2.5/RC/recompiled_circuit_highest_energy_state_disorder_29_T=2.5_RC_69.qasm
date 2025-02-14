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
rz(1.2443378) q[0];
sx q[0];
rz(-0.79255784) q[0];
sx q[0];
rz(-0.25294024) q[0];
rz(-0.77415544) q[1];
sx q[1];
rz(-0.3781265) q[1];
sx q[1];
rz(-2.7546496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16432504) q[0];
sx q[0];
rz(-1.8935793) q[0];
sx q[0];
rz(0.27077814) q[0];
rz(-pi) q[1];
rz(1.4841485) q[2];
sx q[2];
rz(-0.99471015) q[2];
sx q[2];
rz(0.91743166) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36648892) q[1];
sx q[1];
rz(-1.6548043) q[1];
sx q[1];
rz(1.6506399) q[1];
rz(-pi) q[2];
rz(-2.4693842) q[3];
sx q[3];
rz(-2.432193) q[3];
sx q[3];
rz(-1.0649452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0953377) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(-1.7389899) q[2];
rz(1.7272353) q[3];
sx q[3];
rz(-0.71310133) q[3];
sx q[3];
rz(0.49899092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614748) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(2.4312191) q[0];
rz(-1.0164227) q[1];
sx q[1];
rz(-1.2998394) q[1];
sx q[1];
rz(-0.76510915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7920096) q[0];
sx q[0];
rz(-2.5044764) q[0];
sx q[0];
rz(-2.5883915) q[0];
rz(0.76025195) q[2];
sx q[2];
rz(-0.85101247) q[2];
sx q[2];
rz(2.8362897) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5835692) q[1];
sx q[1];
rz(-2.2127732) q[1];
sx q[1];
rz(-0.40972538) q[1];
rz(0.14987544) q[3];
sx q[3];
rz(-1.9016163) q[3];
sx q[3];
rz(0.56043032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1171099) q[2];
sx q[2];
rz(-0.64579248) q[2];
sx q[2];
rz(-2.4369241) q[2];
rz(0.17413983) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(-2.7599938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0524549) q[0];
sx q[0];
rz(-1.313504) q[0];
sx q[0];
rz(-2.254159) q[0];
rz(1.8916091) q[1];
sx q[1];
rz(-1.9307815) q[1];
sx q[1];
rz(1.3608305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3850573) q[0];
sx q[0];
rz(-2.4891315) q[0];
sx q[0];
rz(-2.4908309) q[0];
rz(0.08000441) q[2];
sx q[2];
rz(-1.4334049) q[2];
sx q[2];
rz(-2.4126855) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2584784) q[1];
sx q[1];
rz(-0.7683903) q[1];
sx q[1];
rz(-1.4264631) q[1];
x q[2];
rz(-1.9599846) q[3];
sx q[3];
rz(-1.5312315) q[3];
sx q[3];
rz(0.54007441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7736194) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(1.492929) q[2];
rz(0.44935539) q[3];
sx q[3];
rz(-1.2949233) q[3];
sx q[3];
rz(1.2202643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.087695) q[0];
sx q[0];
rz(-2.5515285) q[0];
sx q[0];
rz(1.4087403) q[0];
rz(2.159481) q[1];
sx q[1];
rz(-1.5928007) q[1];
sx q[1];
rz(1.7353479) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3228673) q[0];
sx q[0];
rz(-1.7850842) q[0];
sx q[0];
rz(-3.0279798) q[0];
rz(2.2135229) q[2];
sx q[2];
rz(-0.065970369) q[2];
sx q[2];
rz(-2.3118491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0413734) q[1];
sx q[1];
rz(-0.58285671) q[1];
sx q[1];
rz(1.1122056) q[1];
rz(-pi) q[2];
rz(2.8729183) q[3];
sx q[3];
rz(-2.2019535) q[3];
sx q[3];
rz(2.1674066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1534319) q[2];
sx q[2];
rz(-1.016523) q[2];
sx q[2];
rz(-2.9856258) q[2];
rz(-2.4912452) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(-3.0273738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0735737) q[0];
sx q[0];
rz(-1.2696215) q[0];
sx q[0];
rz(-0.20075783) q[0];
rz(-1.9643895) q[1];
sx q[1];
rz(-0.8526082) q[1];
sx q[1];
rz(1.9815365) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8070712) q[0];
sx q[0];
rz(-1.40309) q[0];
sx q[0];
rz(1.7981862) q[0];
rz(3.026408) q[2];
sx q[2];
rz(-1.5012001) q[2];
sx q[2];
rz(-1.7382966) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9721891) q[1];
sx q[1];
rz(-2.6693925) q[1];
sx q[1];
rz(0.71581033) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2358642) q[3];
sx q[3];
rz(-2.5865002) q[3];
sx q[3];
rz(2.5344283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3843627) q[2];
sx q[2];
rz(-2.195916) q[2];
sx q[2];
rz(0.45822701) q[2];
rz(0.630817) q[3];
sx q[3];
rz(-1.9061371) q[3];
sx q[3];
rz(0.49921504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.683627) q[0];
sx q[0];
rz(-2.2891335) q[0];
sx q[0];
rz(-0.0068579554) q[0];
rz(-2.9099756) q[1];
sx q[1];
rz(-1.3791142) q[1];
sx q[1];
rz(1.0622271) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4214913) q[0];
sx q[0];
rz(-1.4473697) q[0];
sx q[0];
rz(2.8272259) q[0];
x q[1];
rz(1.7825323) q[2];
sx q[2];
rz(-2.0604302) q[2];
sx q[2];
rz(2.8608866) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1072709) q[1];
sx q[1];
rz(-0.91718972) q[1];
sx q[1];
rz(-1.8922217) q[1];
rz(2.9763664) q[3];
sx q[3];
rz(-0.86841419) q[3];
sx q[3];
rz(-0.14152292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14083938) q[2];
sx q[2];
rz(-1.8893628) q[2];
sx q[2];
rz(1.2678649) q[2];
rz(2.2800692) q[3];
sx q[3];
rz(-2.0778766) q[3];
sx q[3];
rz(-1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93290257) q[0];
sx q[0];
rz(-0.33354315) q[0];
sx q[0];
rz(-2.9260337) q[0];
rz(-2.6453099) q[1];
sx q[1];
rz(-1.9535306) q[1];
sx q[1];
rz(-2.9821679) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2886141) q[0];
sx q[0];
rz(-2.5067177) q[0];
sx q[0];
rz(-1.5663516) q[0];
rz(-pi) q[1];
rz(0.47093289) q[2];
sx q[2];
rz(-1.3976946) q[2];
sx q[2];
rz(1.7488232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9558952) q[1];
sx q[1];
rz(-1.61191) q[1];
sx q[1];
rz(1.3354882) q[1];
rz(-1.0219021) q[3];
sx q[3];
rz(-1.9733784) q[3];
sx q[3];
rz(-2.2233913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7364007) q[2];
sx q[2];
rz(-2.0173732) q[2];
sx q[2];
rz(2.6250725) q[2];
rz(-1.9308331) q[3];
sx q[3];
rz(-1.6885933) q[3];
sx q[3];
rz(1.7621382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6571758) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(2.6013689) q[0];
rz(1.6172488) q[1];
sx q[1];
rz(-2.1397619) q[1];
sx q[1];
rz(-1.8744972) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62477797) q[0];
sx q[0];
rz(-2.5886726) q[0];
sx q[0];
rz(-2.7163158) q[0];
x q[1];
rz(-2.2742496) q[2];
sx q[2];
rz(-2.2144285) q[2];
sx q[2];
rz(-2.1842253) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6798576) q[1];
sx q[1];
rz(-1.9377794) q[1];
sx q[1];
rz(-1.7253897) q[1];
rz(2.2121939) q[3];
sx q[3];
rz(-1.2817973) q[3];
sx q[3];
rz(2.3232587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8871062) q[2];
sx q[2];
rz(-1.1595414) q[2];
sx q[2];
rz(-0.17733388) q[2];
rz(1.0698498) q[3];
sx q[3];
rz(-1.0515352) q[3];
sx q[3];
rz(0.18226084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3607445) q[0];
sx q[0];
rz(-1.1540664) q[0];
sx q[0];
rz(3.0408707) q[0];
rz(-1.1489457) q[1];
sx q[1];
rz(-1.1958242) q[1];
sx q[1];
rz(1.1740059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3927395) q[0];
sx q[0];
rz(-2.1707188) q[0];
sx q[0];
rz(-1.8003182) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4743382) q[2];
sx q[2];
rz(-2.3347179) q[2];
sx q[2];
rz(-2.684456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9575038) q[1];
sx q[1];
rz(-0.042467707) q[1];
sx q[1];
rz(1.0951418) q[1];
rz(0.87971808) q[3];
sx q[3];
rz(-1.4914601) q[3];
sx q[3];
rz(-0.21896958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.164087) q[2];
sx q[2];
rz(-1.6567433) q[2];
sx q[2];
rz(2.738415) q[2];
rz(2.4781503) q[3];
sx q[3];
rz(-2.5622029) q[3];
sx q[3];
rz(-0.0042393953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73096257) q[0];
sx q[0];
rz(-1.0799438) q[0];
sx q[0];
rz(-0.078911111) q[0];
rz(2.2321189) q[1];
sx q[1];
rz(-1.0458922) q[1];
sx q[1];
rz(-0.39631072) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3003259) q[0];
sx q[0];
rz(-1.5864276) q[0];
sx q[0];
rz(-3.1325186) q[0];
x q[1];
rz(1.5780007) q[2];
sx q[2];
rz(-2.3644937) q[2];
sx q[2];
rz(0.33852984) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4034319) q[1];
sx q[1];
rz(-0.53817525) q[1];
sx q[1];
rz(-2.6713283) q[1];
x q[2];
rz(2.2811398) q[3];
sx q[3];
rz(-0.75569587) q[3];
sx q[3];
rz(-0.44482728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4497455) q[2];
sx q[2];
rz(-0.88374603) q[2];
sx q[2];
rz(-0.79908243) q[2];
rz(-0.07829047) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(-2.4962795) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4770724) q[0];
sx q[0];
rz(-1.7417396) q[0];
sx q[0];
rz(-2.59792) q[0];
rz(-1.1935344) q[1];
sx q[1];
rz(-1.2763034) q[1];
sx q[1];
rz(-2.6313849) q[1];
rz(-2.0281744) q[2];
sx q[2];
rz(-1.7241679) q[2];
sx q[2];
rz(0.67571251) q[2];
rz(2.067749) q[3];
sx q[3];
rz(-2.2324149) q[3];
sx q[3];
rz(0.80692337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(6.7232806) q[0];
sx q[0];
rz(6.4203782) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(-0.52991968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13745452) q[0];
sx q[0];
rz(-1.9507427) q[0];
sx q[0];
rz(-3.0259892) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9723188) q[2];
sx q[2];
rz(-1.2520773) q[2];
sx q[2];
rz(-2.4674494) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6119541) q[1];
sx q[1];
rz(-1.0736335) q[1];
sx q[1];
rz(2.0938718) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3210117) q[3];
sx q[3];
rz(-2.3134109) q[3];
sx q[3];
rz(3.0299377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(2.8049862) q[2];
rz(1.6254788) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(1.6158993) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34823725) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(0.021214699) q[0];
rz(1.9477828) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-0.83591998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0854557) q[0];
sx q[0];
rz(-1.5888927) q[0];
sx q[0];
rz(-1.7096814) q[0];
rz(-pi) q[1];
rz(-2.1199273) q[2];
sx q[2];
rz(-1.2163391) q[2];
sx q[2];
rz(0.82565386) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9906209) q[1];
sx q[1];
rz(-0.70578209) q[1];
sx q[1];
rz(1.1433931) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1777982) q[3];
sx q[3];
rz(-2.2713695) q[3];
sx q[3];
rz(0.29004471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8643643) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(1.345984) q[2];
rz(-0.35955444) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(-0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(1.0536449) q[0];
rz(-1.9127649) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-2.704481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489704) q[0];
sx q[0];
rz(-1.4842352) q[0];
sx q[0];
rz(1.2785046) q[0];
rz(2.897981) q[2];
sx q[2];
rz(-1.9677475) q[2];
sx q[2];
rz(1.5543907) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5324085) q[1];
sx q[1];
rz(-0.44645616) q[1];
sx q[1];
rz(-0.32534728) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0760355) q[3];
sx q[3];
rz(-1.6238188) q[3];
sx q[3];
rz(-0.094129063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1221216) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(2.1195228) q[2];
rz(1.2381037) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.5220752) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(3.006382) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-0.19128004) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6387588) q[0];
sx q[0];
rz(-1.6573818) q[0];
sx q[0];
rz(-1.7201352) q[0];
x q[1];
rz(1.26881) q[2];
sx q[2];
rz(-2.3846845) q[2];
sx q[2];
rz(2.7243171) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3595708) q[1];
sx q[1];
rz(-1.4683717) q[1];
sx q[1];
rz(-2.3624079) q[1];
x q[2];
rz(-2.6077765) q[3];
sx q[3];
rz(-1.8752408) q[3];
sx q[3];
rz(2.387405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(-2.130924) q[2];
rz(0.7615532) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(2.9060569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33870944) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(-0.55661911) q[0];
rz(0.11511766) q[1];
sx q[1];
rz(-1.3373673) q[1];
sx q[1];
rz(0.97250485) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9113377) q[0];
sx q[0];
rz(-1.6382434) q[0];
sx q[0];
rz(-3.0393242) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82053484) q[2];
sx q[2];
rz(-1.8168601) q[2];
sx q[2];
rz(2.8503502) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3481808) q[1];
sx q[1];
rz(-0.23985292) q[1];
sx q[1];
rz(2.1972149) q[1];
x q[2];
rz(-3.1049018) q[3];
sx q[3];
rz(-0.9551691) q[3];
sx q[3];
rz(2.7021367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82289034) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.5931607) q[2];
rz(-1.3657773) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(0.95388609) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71972972) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(-1.7156037) q[0];
rz(2.0772207) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(-0.37429601) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4842589) q[0];
sx q[0];
rz(-1.9801635) q[0];
sx q[0];
rz(1.9166458) q[0];
x q[1];
rz(-2.6815368) q[2];
sx q[2];
rz(-2.6003777) q[2];
sx q[2];
rz(-0.40749007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.018505521) q[1];
sx q[1];
rz(-1.0235041) q[1];
sx q[1];
rz(-1.2612543) q[1];
x q[2];
rz(-2.3338942) q[3];
sx q[3];
rz(-0.85113111) q[3];
sx q[3];
rz(-0.52586517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2465683) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(2.1833615) q[2];
rz(-2.9124177) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-0.36528698) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(0.90674415) q[0];
rz(-2.0523741) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.3100756) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.037127) q[0];
sx q[0];
rz(-1.6220399) q[0];
sx q[0];
rz(-0.38802223) q[0];
x q[1];
rz(-1.9016978) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(-2.9943525) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.29855) q[1];
sx q[1];
rz(-1.3805461) q[1];
sx q[1];
rz(-0.33903867) q[1];
rz(-pi) q[2];
rz(-0.044192627) q[3];
sx q[3];
rz(-0.5558388) q[3];
sx q[3];
rz(2.1293872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.92581302) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(1.6581992) q[2];
rz(0.27967927) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(-0.057597615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9782372) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(-2.9220007) q[0];
rz(2.638468) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(0.84987744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69707623) q[0];
sx q[0];
rz(-2.8651617) q[0];
sx q[0];
rz(0.96093775) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88768994) q[2];
sx q[2];
rz(-1.6437093) q[2];
sx q[2];
rz(2.2788252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0975768) q[1];
sx q[1];
rz(-0.8102639) q[1];
sx q[1];
rz(3.0939328) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1342808) q[3];
sx q[3];
rz(-0.28537649) q[3];
sx q[3];
rz(2.7989822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5510817) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(1.760651) q[2];
rz(0.75602174) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(-0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6324156) q[0];
sx q[0];
rz(-0.89288765) q[0];
sx q[0];
rz(2.7365622) q[0];
rz(0.45267725) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(-1.8639494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4687846) q[0];
sx q[0];
rz(-1.0793669) q[0];
sx q[0];
rz(-1.3093033) q[0];
rz(-pi) q[1];
rz(-1.0870886) q[2];
sx q[2];
rz(-0.93855575) q[2];
sx q[2];
rz(1.5916923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7840522) q[1];
sx q[1];
rz(-1.6077542) q[1];
sx q[1];
rz(1.1251015) q[1];
rz(-pi) q[2];
rz(2.7550335) q[3];
sx q[3];
rz(-2.2694526) q[3];
sx q[3];
rz(0.78702918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3383125) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(-2.7837616) q[2];
rz(1.7221649) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2232067) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(3.0293368) q[0];
rz(2.2414801) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(0.21044883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8716547) q[0];
sx q[0];
rz(-2.6007814) q[0];
sx q[0];
rz(-2.7842115) q[0];
x q[1];
rz(1.9333282) q[2];
sx q[2];
rz(-2.7047485) q[2];
sx q[2];
rz(0.7561965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3887742) q[1];
sx q[1];
rz(-0.37467271) q[1];
sx q[1];
rz(-2.545536) q[1];
x q[2];
rz(0.25037346) q[3];
sx q[3];
rz(-1.3406521) q[3];
sx q[3];
rz(1.2311414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18008733) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(-1.2735584) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(-2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56959854) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(2.3251484) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(-2.0408761) q[2];
sx q[2];
rz(-1.7190949) q[2];
sx q[2];
rz(0.20863056) q[2];
rz(2.3002426) q[3];
sx q[3];
rz(-1.8302866) q[3];
sx q[3];
rz(2.4091099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

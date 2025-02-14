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
rz(0.1641195) q[0];
sx q[0];
rz(-2.1174705) q[0];
sx q[0];
rz(0.48409387) q[0];
rz(-0.76264277) q[1];
sx q[1];
rz(-1.5781382) q[1];
sx q[1];
rz(0.054952316) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2919827) q[0];
sx q[0];
rz(-1.6020444) q[0];
sx q[0];
rz(0.90051767) q[0];
rz(1.820716) q[2];
sx q[2];
rz(-1.6276161) q[2];
sx q[2];
rz(-2.394258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72803264) q[1];
sx q[1];
rz(-1.7031341) q[1];
sx q[1];
rz(1.669426) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1922429) q[3];
sx q[3];
rz(-1.8481405) q[3];
sx q[3];
rz(-0.30888939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83076465) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(0.27466276) q[2];
rz(0.87485391) q[3];
sx q[3];
rz(-1.091205) q[3];
sx q[3];
rz(0.27354512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
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
rz(2.4817568) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(-0.39875317) q[0];
rz(-2.8975471) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(-1.1057378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53351346) q[0];
sx q[0];
rz(-1.7987804) q[0];
sx q[0];
rz(-0.45324676) q[0];
rz(-pi) q[1];
rz(-2.2244338) q[2];
sx q[2];
rz(-1.0984525) q[2];
sx q[2];
rz(-1.2231959) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2011021) q[1];
sx q[1];
rz(-0.0041714287) q[1];
sx q[1];
rz(2.4524686) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12064528) q[3];
sx q[3];
rz(-1.835726) q[3];
sx q[3];
rz(2.210833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70046052) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(2.0665118) q[2];
rz(2.4454146) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.8657846) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(2.5221141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86176819) q[0];
sx q[0];
rz(-0.048088308) q[0];
sx q[0];
rz(-2.2347941) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9335494) q[2];
sx q[2];
rz(-3.0100076) q[2];
sx q[2];
rz(-3.0818444) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3302096) q[1];
sx q[1];
rz(-2.7683677) q[1];
sx q[1];
rz(0.95287816) q[1];
x q[2];
rz(-2.3893395) q[3];
sx q[3];
rz(-0.76365389) q[3];
sx q[3];
rz(-1.4925166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2565903) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(-3.0510862) q[2];
rz(-2.6835486) q[3];
sx q[3];
rz(-2.6543255) q[3];
sx q[3];
rz(-1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8103545) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(-0.93165398) q[0];
rz(-2.9725507) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(1.0394675) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19379481) q[0];
sx q[0];
rz(-0.76776869) q[0];
sx q[0];
rz(2.7916273) q[0];
rz(-pi) q[1];
rz(1.3170856) q[2];
sx q[2];
rz(-1.3519962) q[2];
sx q[2];
rz(-1.7148866) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89371496) q[1];
sx q[1];
rz(-1.3637241) q[1];
sx q[1];
rz(-0.061092579) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28431953) q[3];
sx q[3];
rz(-1.7756724) q[3];
sx q[3];
rz(-2.7903872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6197551) q[2];
sx q[2];
rz(-2.0401185) q[2];
sx q[2];
rz(-2.6452765) q[2];
rz(0.79658341) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(1.0930141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26688823) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(0.97187483) q[0];
rz(3.019849) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(-1.3294539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0861391) q[0];
sx q[0];
rz(-0.3124961) q[0];
sx q[0];
rz(1.7715095) q[0];
x q[1];
rz(2.0636286) q[2];
sx q[2];
rz(-1.9415641) q[2];
sx q[2];
rz(-0.28647067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7763627) q[1];
sx q[1];
rz(-0.91971469) q[1];
sx q[1];
rz(1.0427703) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6498783) q[3];
sx q[3];
rz(-1.5800161) q[3];
sx q[3];
rz(0.10445933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10151265) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(-2.9902048) q[2];
rz(2.3769412) q[3];
sx q[3];
rz(-2.305856) q[3];
sx q[3];
rz(-1.5966655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695892) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(-2.0701011) q[0];
rz(-2.6103861) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(-2.5154617) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6036442) q[0];
sx q[0];
rz(-1.5755565) q[0];
sx q[0];
rz(-1.5877923) q[0];
x q[1];
rz(-1.3254204) q[2];
sx q[2];
rz(-0.98773709) q[2];
sx q[2];
rz(2.8193605) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31589139) q[1];
sx q[1];
rz(-2.4448312) q[1];
sx q[1];
rz(2.6543762) q[1];
rz(0.71685426) q[3];
sx q[3];
rz(-0.62590137) q[3];
sx q[3];
rz(-0.38960534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7834187) q[2];
sx q[2];
rz(-2.5219315) q[2];
sx q[2];
rz(1.0142856) q[2];
rz(-1.8861534) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(2.0881418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805873) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(-2.2210806) q[0];
rz(-0.11407425) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(2.0106409) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81984542) q[0];
sx q[0];
rz(-1.5764569) q[0];
sx q[0];
rz(-0.59771718) q[0];
rz(-0.80318309) q[2];
sx q[2];
rz(-0.72942643) q[2];
sx q[2];
rz(2.8620697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8495226) q[1];
sx q[1];
rz(-0.91150586) q[1];
sx q[1];
rz(-0.16772049) q[1];
rz(-pi) q[2];
rz(-3.00131) q[3];
sx q[3];
rz(-0.52061235) q[3];
sx q[3];
rz(-1.8050107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(0.40880173) q[2];
rz(0.18320228) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278397) q[0];
sx q[0];
rz(-0.04700679) q[0];
sx q[0];
rz(2.8507932) q[0];
rz(0.94789061) q[1];
sx q[1];
rz(-2.6746076) q[1];
sx q[1];
rz(1.1999757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5644218) q[0];
sx q[0];
rz(-1.0438215) q[0];
sx q[0];
rz(2.0354664) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84917111) q[2];
sx q[2];
rz(-0.90540041) q[2];
sx q[2];
rz(-0.1553436) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6281451) q[1];
sx q[1];
rz(-1.8643059) q[1];
sx q[1];
rz(-0.048444466) q[1];
rz(3.0407716) q[3];
sx q[3];
rz(-1.7530879) q[3];
sx q[3];
rz(1.3829447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7802508) q[2];
sx q[2];
rz(-1.6929408) q[2];
sx q[2];
rz(1.7928436) q[2];
rz(-1.5032984) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(-2.9564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0272442) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(1.0741023) q[0];
rz(-2.7117924) q[1];
sx q[1];
rz(-1.0440412) q[1];
sx q[1];
rz(2.6780186) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0811679) q[0];
sx q[0];
rz(-1.3827818) q[0];
sx q[0];
rz(3.0449165) q[0];
x q[1];
rz(-2.8797382) q[2];
sx q[2];
rz(-1.3758278) q[2];
sx q[2];
rz(-1.5972114) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.025109865) q[1];
sx q[1];
rz(-1.1949202) q[1];
sx q[1];
rz(0.89565887) q[1];
x q[2];
rz(-1.4809247) q[3];
sx q[3];
rz(-1.0139795) q[3];
sx q[3];
rz(-1.1695605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27434906) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(-2.3331433) q[2];
rz(2.6570053) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(0.43420473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3633858) q[0];
sx q[0];
rz(-0.45553842) q[0];
sx q[0];
rz(-2.0090012) q[0];
rz(-3.1032108) q[1];
sx q[1];
rz(-1.8033586) q[1];
sx q[1];
rz(-2.8448232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0050186) q[0];
sx q[0];
rz(-1.8474192) q[0];
sx q[0];
rz(2.121595) q[0];
rz(-pi) q[1];
rz(0.77905853) q[2];
sx q[2];
rz(-0.59528661) q[2];
sx q[2];
rz(0.99814683) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0368288) q[1];
sx q[1];
rz(-1.4955758) q[1];
sx q[1];
rz(0.22237088) q[1];
x q[2];
rz(1.1007916) q[3];
sx q[3];
rz(-0.74277011) q[3];
sx q[3];
rz(2.7239885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.177629) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(0.56619823) q[2];
rz(-1.1244134) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(-1.9521149) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2531256) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(2.879907) q[1];
sx q[1];
rz(-1.9410004) q[1];
sx q[1];
rz(0.95071361) q[1];
rz(-0.08200866) q[2];
sx q[2];
rz(-1.5951984) q[2];
sx q[2];
rz(-2.3328608) q[2];
rz(2.2903743) q[3];
sx q[3];
rz(-0.85243445) q[3];
sx q[3];
rz(1.7029521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(2.6296122) q[0];
sx q[0];
rz(-0.57096243) q[0];
sx q[0];
rz(2.9538739) q[0];
rz(-2.3277148) q[1];
sx q[1];
rz(-1.4717646) q[1];
sx q[1];
rz(1.8522813) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16351249) q[0];
sx q[0];
rz(-0.31766444) q[0];
sx q[0];
rz(-0.72084041) q[0];
rz(-pi) q[1];
rz(-3.0953636) q[2];
sx q[2];
rz(-2.5704489) q[2];
sx q[2];
rz(0.91396224) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.493452) q[1];
sx q[1];
rz(-1.3228184) q[1];
sx q[1];
rz(2.6676755) q[1];
rz(-2.1407152) q[3];
sx q[3];
rz(-2.7288247) q[3];
sx q[3];
rz(1.643553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2671555) q[2];
sx q[2];
rz(-2.0384553) q[2];
sx q[2];
rz(2.3647986) q[2];
rz(-1.3022425) q[3];
sx q[3];
rz(-1.6603419) q[3];
sx q[3];
rz(-1.9384025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60748196) q[0];
sx q[0];
rz(-0.48180875) q[0];
sx q[0];
rz(2.1433461) q[0];
rz(-2.7718995) q[1];
sx q[1];
rz(-1.0414711) q[1];
sx q[1];
rz(2.0236156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6551483) q[0];
sx q[0];
rz(-1.0308497) q[0];
sx q[0];
rz(2.6686431) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68685617) q[2];
sx q[2];
rz(-2.0389028) q[2];
sx q[2];
rz(1.8091701) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79827944) q[1];
sx q[1];
rz(-1.717853) q[1];
sx q[1];
rz(-1.1215854) q[1];
rz(-pi) q[2];
rz(0.93922521) q[3];
sx q[3];
rz(-1.7895921) q[3];
sx q[3];
rz(-2.95524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.86576858) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(-2.7562874) q[2];
rz(-0.038711874) q[3];
sx q[3];
rz(-2.7888515) q[3];
sx q[3];
rz(2.501131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63967079) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(-1.3013526) q[0];
rz(3.0838857) q[1];
sx q[1];
rz(-0.4464018) q[1];
sx q[1];
rz(-1.3267964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2952174) q[0];
sx q[0];
rz(-2.1070988) q[0];
sx q[0];
rz(0.55732507) q[0];
rz(-pi) q[1];
rz(-2.9811601) q[2];
sx q[2];
rz(-1.7084873) q[2];
sx q[2];
rz(1.4169803) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2542299) q[1];
sx q[1];
rz(-2.0944203) q[1];
sx q[1];
rz(1.2382522) q[1];
rz(-pi) q[2];
rz(2.7884198) q[3];
sx q[3];
rz(-2.2614217) q[3];
sx q[3];
rz(-1.4598802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79746276) q[2];
sx q[2];
rz(-2.3064488) q[2];
sx q[2];
rz(2.3878035) q[2];
rz(-1.7720743) q[3];
sx q[3];
rz(-1.6794208) q[3];
sx q[3];
rz(2.4863825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.9264483) q[0];
sx q[0];
rz(-1.0855874) q[0];
sx q[0];
rz(1.7468859) q[0];
rz(1.9649547) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(0.36697695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4876539) q[0];
sx q[0];
rz(-1.5040845) q[0];
sx q[0];
rz(-3.0848461) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4796333) q[2];
sx q[2];
rz(-0.41482224) q[2];
sx q[2];
rz(-2.7560459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9865468) q[1];
sx q[1];
rz(-1.9268225) q[1];
sx q[1];
rz(-1.3638391) q[1];
rz(1.4848726) q[3];
sx q[3];
rz(-1.9367227) q[3];
sx q[3];
rz(0.56473063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7299812) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(-1.270594) q[2];
rz(-0.96108428) q[3];
sx q[3];
rz(-1.7226487) q[3];
sx q[3];
rz(-3.0911176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5807895) q[0];
sx q[0];
rz(-0.92629782) q[0];
sx q[0];
rz(1.3128989) q[0];
rz(-1.987223) q[1];
sx q[1];
rz(-1.9998974) q[1];
sx q[1];
rz(-0.55783522) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3274021) q[0];
sx q[0];
rz(-2.0923776) q[0];
sx q[0];
rz(-0.5970229) q[0];
x q[1];
rz(7/(5*pi)) q[2];
sx q[2];
rz(-0.17765309) q[2];
sx q[2];
rz(-0.38933094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3250998) q[1];
sx q[1];
rz(-1.268728) q[1];
sx q[1];
rz(2.3343071) q[1];
rz(-0.29739012) q[3];
sx q[3];
rz(-0.47978401) q[3];
sx q[3];
rz(-0.23272091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4981726) q[2];
sx q[2];
rz(-1.2926481) q[2];
sx q[2];
rz(-0.84929973) q[2];
rz(0.15527209) q[3];
sx q[3];
rz(-0.97584358) q[3];
sx q[3];
rz(0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3575386) q[0];
sx q[0];
rz(-2.5548866) q[0];
sx q[0];
rz(-0.46911711) q[0];
rz(-0.47438374) q[1];
sx q[1];
rz(-2.3553039) q[1];
sx q[1];
rz(-2.3862086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84635272) q[0];
sx q[0];
rz(-1.8420514) q[0];
sx q[0];
rz(-0.004581906) q[0];
rz(1.239849) q[2];
sx q[2];
rz(-2.1734383) q[2];
sx q[2];
rz(-1.861426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0363706) q[1];
sx q[1];
rz(-1.7835938) q[1];
sx q[1];
rz(0.12963055) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83306082) q[3];
sx q[3];
rz(-1.2445881) q[3];
sx q[3];
rz(1.2556094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0581806) q[2];
sx q[2];
rz(-1.3231134) q[2];
sx q[2];
rz(-2.9537436) q[2];
rz(1.2457054) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(2.0511621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2096527) q[0];
sx q[0];
rz(-1.8745475) q[0];
sx q[0];
rz(2.7506822) q[0];
rz(-0.93005013) q[1];
sx q[1];
rz(-1.4314194) q[1];
sx q[1];
rz(-1.3040868) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5374889) q[0];
sx q[0];
rz(-1.2199645) q[0];
sx q[0];
rz(3.119191) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9043998) q[2];
sx q[2];
rz(-2.790934) q[2];
sx q[2];
rz(1.2192977) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89829373) q[1];
sx q[1];
rz(-2.0212681) q[1];
sx q[1];
rz(2.9551278) q[1];
rz(-0.061873925) q[3];
sx q[3];
rz(-2.8447667) q[3];
sx q[3];
rz(-0.59943953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8730674) q[2];
sx q[2];
rz(-2.8748942) q[2];
sx q[2];
rz(3.0282057) q[2];
rz(3.125627) q[3];
sx q[3];
rz(-1.6458052) q[3];
sx q[3];
rz(-2.9769843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78405821) q[0];
sx q[0];
rz(-1.6679732) q[0];
sx q[0];
rz(0.12406021) q[0];
rz(-0.1217753) q[1];
sx q[1];
rz(-0.75141326) q[1];
sx q[1];
rz(1.442499) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6704491) q[0];
sx q[0];
rz(-1.6713977) q[0];
sx q[0];
rz(1.9557682) q[0];
x q[1];
rz(2.997918) q[2];
sx q[2];
rz(-2.8468003) q[2];
sx q[2];
rz(-0.93570342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9297819) q[1];
sx q[1];
rz(-2.7719797) q[1];
sx q[1];
rz(-1.1951642) q[1];
rz(-2.8644531) q[3];
sx q[3];
rz(-2.3921514) q[3];
sx q[3];
rz(-2.998874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1759935) q[2];
sx q[2];
rz(-1.7815353) q[2];
sx q[2];
rz(2.3967801) q[2];
rz(-2.2143769) q[3];
sx q[3];
rz(-1.6330481) q[3];
sx q[3];
rz(-2.0492699) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6532779) q[0];
sx q[0];
rz(-1.7019685) q[0];
sx q[0];
rz(-0.22536817) q[0];
rz(-1.1124181) q[1];
sx q[1];
rz(-2.367159) q[1];
sx q[1];
rz(-2.2427799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5170636) q[0];
sx q[0];
rz(-1.0072021) q[0];
sx q[0];
rz(-0.64080142) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5139047) q[2];
sx q[2];
rz(-1.9348839) q[2];
sx q[2];
rz(2.4533426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73279335) q[1];
sx q[1];
rz(-2.2299754) q[1];
sx q[1];
rz(-2.885347) q[1];
rz(-pi) q[2];
x q[2];
rz(3.131116) q[3];
sx q[3];
rz(-1.8636522) q[3];
sx q[3];
rz(1.7550857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0972458) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(-2.782235) q[2];
rz(0.25092956) q[3];
sx q[3];
rz(-1.072262) q[3];
sx q[3];
rz(2.3738764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3971685) q[0];
sx q[0];
rz(-0.1305307) q[0];
sx q[0];
rz(-0.8771483) q[0];
rz(1.5646704) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(-2.3559779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1164074) q[0];
sx q[0];
rz(-1.8078363) q[0];
sx q[0];
rz(-2.9171506) q[0];
x q[1];
rz(-1.142611) q[2];
sx q[2];
rz(-2.0642363) q[2];
sx q[2];
rz(2.7827415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4686376) q[1];
sx q[1];
rz(-1.6412927) q[1];
sx q[1];
rz(-2.7967909) q[1];
rz(-pi) q[2];
rz(1.7212058) q[3];
sx q[3];
rz(-2.0724943) q[3];
sx q[3];
rz(0.13817638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-2.3515297) q[2];
sx q[2];
rz(2.4033974) q[2];
rz(-0.92292845) q[3];
sx q[3];
rz(-1.3083369) q[3];
sx q[3];
rz(-0.2963399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34687635) q[0];
sx q[0];
rz(-1.7807757) q[0];
sx q[0];
rz(-2.1697252) q[0];
rz(2.234266) q[1];
sx q[1];
rz(-2.0569888) q[1];
sx q[1];
rz(2.8203698) q[1];
rz(2.7493334) q[2];
sx q[2];
rz(-1.0259368) q[2];
sx q[2];
rz(-2.1403014) q[2];
rz(-0.76413705) q[3];
sx q[3];
rz(-2.4485179) q[3];
sx q[3];
rz(-2.5701523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-2.0353844) q[0];
sx q[0];
rz(2.455403) q[0];
sx q[0];
rz(10.578293) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(-2.6931813) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9759244) q[0];
sx q[0];
rz(-2.2184555) q[0];
sx q[0];
rz(1.5041385) q[0];
rz(-pi) q[1];
rz(1.6719867) q[2];
sx q[2];
rz(-1.8253583) q[2];
sx q[2];
rz(2.2806666) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37473512) q[1];
sx q[1];
rz(-0.85588662) q[1];
sx q[1];
rz(-0.84918569) q[1];
rz(-pi) q[2];
rz(1.6187864) q[3];
sx q[3];
rz(-2.0715044) q[3];
sx q[3];
rz(-1.9013311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6260234) q[2];
sx q[2];
rz(-1.7531351) q[2];
sx q[2];
rz(2.5173729) q[2];
rz(0.37912399) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(-0.24851255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17495951) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(2.1061184) q[0];
rz(1.2738719) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(-0.48286352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0014526) q[0];
sx q[0];
rz(-1.4743544) q[0];
sx q[0];
rz(-1.6167377) q[0];
x q[1];
rz(2.9752983) q[2];
sx q[2];
rz(-1.6876843) q[2];
sx q[2];
rz(0.4695732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.711593) q[1];
sx q[1];
rz(-1.5801799) q[1];
sx q[1];
rz(-1.1389772) q[1];
rz(-pi) q[2];
rz(-1.6522406) q[3];
sx q[3];
rz(-1.2508498) q[3];
sx q[3];
rz(0.68405747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.024293385) q[2];
sx q[2];
rz(-1.5017941) q[2];
sx q[2];
rz(-1.8710322) q[2];
rz(-0.67000669) q[3];
sx q[3];
rz(-0.62205258) q[3];
sx q[3];
rz(0.89699927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73827493) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(2.4025412) q[0];
rz(-2.1801379) q[1];
sx q[1];
rz(-0.32197222) q[1];
sx q[1];
rz(-2.3846073) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8563499) q[0];
sx q[0];
rz(-2.3654571) q[0];
sx q[0];
rz(1.9479284) q[0];
x q[1];
rz(-0.44481014) q[2];
sx q[2];
rz(-2.399657) q[2];
sx q[2];
rz(0.43659376) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18338246) q[1];
sx q[1];
rz(-0.64309769) q[1];
sx q[1];
rz(-0.22497589) q[1];
rz(1.3540165) q[3];
sx q[3];
rz(-0.31459537) q[3];
sx q[3];
rz(-1.5217239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5028533) q[2];
sx q[2];
rz(-0.91705489) q[2];
sx q[2];
rz(0.70510954) q[2];
rz(2.8201568) q[3];
sx q[3];
rz(-2.1240081) q[3];
sx q[3];
rz(-0.41079918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(3.131677) q[0];
sx q[0];
rz(-0.83808815) q[0];
sx q[0];
rz(1.9158844) q[0];
rz(-1.2340087) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(-0.066224901) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4915337) q[0];
sx q[0];
rz(-0.26187944) q[0];
sx q[0];
rz(-1.9240379) q[0];
rz(-pi) q[1];
rz(-0.80133665) q[2];
sx q[2];
rz(-1.3532012) q[2];
sx q[2];
rz(0.8720397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60053101) q[1];
sx q[1];
rz(-2.5503073) q[1];
sx q[1];
rz(0.97452428) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6965167) q[3];
sx q[3];
rz(-1.9341521) q[3];
sx q[3];
rz(-1.1896127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0951198) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(2.7009916) q[2];
rz(2.1353841) q[3];
sx q[3];
rz(-1.7677842) q[3];
sx q[3];
rz(0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7202268) q[0];
sx q[0];
rz(-0.81291968) q[0];
sx q[0];
rz(-1.6356069) q[0];
rz(1.6141363) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(1.6857326) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57350006) q[0];
sx q[0];
rz(-1.4254036) q[0];
sx q[0];
rz(0.83515867) q[0];
rz(-1.6873932) q[2];
sx q[2];
rz(-1.3449838) q[2];
sx q[2];
rz(1.1979529) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.887874) q[1];
sx q[1];
rz(-1.591133) q[1];
sx q[1];
rz(-0.65171297) q[1];
rz(-pi) q[2];
rz(0.21946986) q[3];
sx q[3];
rz(-0.61958909) q[3];
sx q[3];
rz(2.0989024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8684034) q[2];
sx q[2];
rz(-1.5634147) q[2];
sx q[2];
rz(0.91147649) q[2];
rz(2.2677926) q[3];
sx q[3];
rz(-2.3995212) q[3];
sx q[3];
rz(-3.1349283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28354302) q[0];
sx q[0];
rz(-0.25074211) q[0];
sx q[0];
rz(0.13033303) q[0];
rz(2.6538972) q[1];
sx q[1];
rz(-0.54134381) q[1];
sx q[1];
rz(-2.3977051) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.795563) q[0];
sx q[0];
rz(-1.6206121) q[0];
sx q[0];
rz(1.4883785) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7623474) q[2];
sx q[2];
rz(-2.5073176) q[2];
sx q[2];
rz(-2.3630362) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53092693) q[1];
sx q[1];
rz(-2.1082889) q[1];
sx q[1];
rz(-2.8100138) q[1];
x q[2];
rz(-1.3906995) q[3];
sx q[3];
rz(-1.0922179) q[3];
sx q[3];
rz(0.046600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.1193715) q[2];
sx q[2];
rz(-2.2890942) q[2];
sx q[2];
rz(0.96088299) q[2];
rz(2.7458701) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(-3.0254288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845067) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(1.047026) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.8641169) q[1];
sx q[1];
rz(-0.39271694) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0020802845) q[0];
sx q[0];
rz(-0.89186984) q[0];
sx q[0];
rz(-0.2261136) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95177841) q[2];
sx q[2];
rz(-1.5238395) q[2];
sx q[2];
rz(2.6359735) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8046569) q[1];
sx q[1];
rz(-1.5976904) q[1];
sx q[1];
rz(-2.8662445) q[1];
rz(-2.8128871) q[3];
sx q[3];
rz(-1.3188601) q[3];
sx q[3];
rz(-3.1303034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0543694) q[2];
sx q[2];
rz(-1.9471709) q[2];
sx q[2];
rz(1.8325904) q[2];
rz(1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(-0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8488309) q[0];
sx q[0];
rz(-1.6852385) q[0];
sx q[0];
rz(0.30717474) q[0];
rz(-1.81709) q[1];
sx q[1];
rz(-0.38506404) q[1];
sx q[1];
rz(-0.90075341) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9884315) q[0];
sx q[0];
rz(-1.5539196) q[0];
sx q[0];
rz(2.8794226) q[0];
rz(-pi) q[1];
rz(1.0020761) q[2];
sx q[2];
rz(-2.4042712) q[2];
sx q[2];
rz(1.7676644) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5430205) q[1];
sx q[1];
rz(-1.9024339) q[1];
sx q[1];
rz(0.33473067) q[1];
rz(-pi) q[2];
rz(-1.4291271) q[3];
sx q[3];
rz(-0.9441388) q[3];
sx q[3];
rz(3.0350445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8073392) q[2];
sx q[2];
rz(-3.0901577) q[2];
sx q[2];
rz(-2.5267498) q[2];
rz(-1.0963415) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(-0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39356247) q[0];
sx q[0];
rz(-0.63069558) q[0];
sx q[0];
rz(-0.50931859) q[0];
rz(-0.18109426) q[1];
sx q[1];
rz(-1.7362005) q[1];
sx q[1];
rz(0.96955713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5405131) q[0];
sx q[0];
rz(-0.76427312) q[0];
sx q[0];
rz(0.22354062) q[0];
x q[1];
rz(-1.4793899) q[2];
sx q[2];
rz(-0.50894032) q[2];
sx q[2];
rz(1.3846874) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73652112) q[1];
sx q[1];
rz(-2.7847291) q[1];
sx q[1];
rz(-0.31330152) q[1];
rz(0.84253715) q[3];
sx q[3];
rz(-1.8368097) q[3];
sx q[3];
rz(-2.0913948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2481044) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-2.060037) q[2];
rz(0.23165101) q[3];
sx q[3];
rz(-1.3737498) q[3];
sx q[3];
rz(-2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.76081) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(0.64714062) q[0];
rz(-0.45516792) q[1];
sx q[1];
rz(-1.7166694) q[1];
sx q[1];
rz(-0.761935) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9058911) q[0];
sx q[0];
rz(-2.3755223) q[0];
sx q[0];
rz(0.42278843) q[0];
x q[1];
rz(2.2438887) q[2];
sx q[2];
rz(-0.87397777) q[2];
sx q[2];
rz(-0.76329939) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6651407) q[1];
sx q[1];
rz(-2.2287205) q[1];
sx q[1];
rz(1.5275147) q[1];
x q[2];
rz(-2.0382138) q[3];
sx q[3];
rz(-1.5639135) q[3];
sx q[3];
rz(2.9228743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0882988) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(2.2376412) q[2];
rz(-1.6186742) q[3];
sx q[3];
rz(-1.0205597) q[3];
sx q[3];
rz(-1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.224613) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(-0.12679535) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(-1.1197208) q[2];
sx q[2];
rz(-1.3430165) q[2];
sx q[2];
rz(-2.7271885) q[2];
rz(-1.5736754) q[3];
sx q[3];
rz(-1.1408014) q[3];
sx q[3];
rz(-0.037527966) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

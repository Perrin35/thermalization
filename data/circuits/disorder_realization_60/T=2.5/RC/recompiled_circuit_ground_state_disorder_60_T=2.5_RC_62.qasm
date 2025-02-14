OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.056501) q[0];
sx q[0];
rz(-2.8488475) q[0];
sx q[0];
rz(-2.2226287) q[0];
rz(-0.38209823) q[1];
sx q[1];
rz(-2.9736019) q[1];
sx q[1];
rz(1.1408495) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1699693) q[0];
sx q[0];
rz(-2.9282018) q[0];
sx q[0];
rz(2.1687228) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2873285) q[2];
sx q[2];
rz(-1.8745443) q[2];
sx q[2];
rz(-1.8423353) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6534868) q[1];
sx q[1];
rz(-1.7035023) q[1];
sx q[1];
rz(2.7567171) q[1];
x q[2];
rz(3.0039722) q[3];
sx q[3];
rz(-1.1351661) q[3];
sx q[3];
rz(2.9916258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69316489) q[2];
sx q[2];
rz(-2.0824671) q[2];
sx q[2];
rz(1.1761752) q[2];
rz(-0.042393427) q[3];
sx q[3];
rz(-1.3375125) q[3];
sx q[3];
rz(-0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6723044) q[0];
sx q[0];
rz(-2.7870218) q[0];
sx q[0];
rz(0.18445036) q[0];
rz(-0.09672673) q[1];
sx q[1];
rz(-1.6177142) q[1];
sx q[1];
rz(1.2299889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0320659) q[0];
sx q[0];
rz(-1.3680887) q[0];
sx q[0];
rz(-0.37287449) q[0];
x q[1];
rz(0.076893496) q[2];
sx q[2];
rz(-2.580759) q[2];
sx q[2];
rz(0.70181134) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7736771) q[1];
sx q[1];
rz(-1.333548) q[1];
sx q[1];
rz(1.3648454) q[1];
rz(-pi) q[2];
rz(3.1026479) q[3];
sx q[3];
rz(-1.9488261) q[3];
sx q[3];
rz(1.259144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7868598) q[2];
sx q[2];
rz(-2.730098) q[2];
sx q[2];
rz(-3.1301609) q[2];
rz(1.8276021) q[3];
sx q[3];
rz(-1.3649789) q[3];
sx q[3];
rz(-2.438681) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74362022) q[0];
sx q[0];
rz(-0.13298661) q[0];
sx q[0];
rz(-2.5308894) q[0];
rz(-0.01344219) q[1];
sx q[1];
rz(-1.2862658) q[1];
sx q[1];
rz(0.59649831) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4152348) q[0];
sx q[0];
rz(-1.4044552) q[0];
sx q[0];
rz(-2.1852058) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99278583) q[2];
sx q[2];
rz(-1.3137378) q[2];
sx q[2];
rz(1.7456733) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92745078) q[1];
sx q[1];
rz(-0.79793707) q[1];
sx q[1];
rz(0.44254704) q[1];
x q[2];
rz(-0.032993837) q[3];
sx q[3];
rz(-1.0171431) q[3];
sx q[3];
rz(-0.43137303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67768031) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(0.00027351969) q[2];
rz(1.6655946) q[3];
sx q[3];
rz(-0.6690343) q[3];
sx q[3];
rz(-3.0345548) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6855455) q[0];
sx q[0];
rz(-2.9750329) q[0];
sx q[0];
rz(-1.9786932) q[0];
rz(-0.97745013) q[1];
sx q[1];
rz(-1.5403427) q[1];
sx q[1];
rz(1.5553442) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9396203) q[0];
sx q[0];
rz(-3.1292136) q[0];
sx q[0];
rz(-2.2511961) q[0];
rz(-0.78537099) q[2];
sx q[2];
rz(-1.5391239) q[2];
sx q[2];
rz(-2.9065064) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0229637) q[1];
sx q[1];
rz(-1.5852155) q[1];
sx q[1];
rz(-0.49683907) q[1];
x q[2];
rz(2.8951696) q[3];
sx q[3];
rz(-1.0895673) q[3];
sx q[3];
rz(-2.6068166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.9419452) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(-2.9072348) q[2];
rz(2.6394081) q[3];
sx q[3];
rz(-1.0415223) q[3];
sx q[3];
rz(-2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31768826) q[0];
sx q[0];
rz(-2.3264139) q[0];
sx q[0];
rz(-1.748388) q[0];
rz(-0.69270095) q[1];
sx q[1];
rz(-2.0557978) q[1];
sx q[1];
rz(1.0903953) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6818974) q[0];
sx q[0];
rz(-0.57006627) q[0];
sx q[0];
rz(1.1718114) q[0];
rz(-2.1887921) q[2];
sx q[2];
rz(-2.0119466) q[2];
sx q[2];
rz(-1.49681) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.322256) q[1];
sx q[1];
rz(-0.44274477) q[1];
sx q[1];
rz(1.2135452) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0308068) q[3];
sx q[3];
rz(-1.4362236) q[3];
sx q[3];
rz(2.7887044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56110567) q[2];
sx q[2];
rz(-1.9876391) q[2];
sx q[2];
rz(-0.53606501) q[2];
rz(2.8481893) q[3];
sx q[3];
rz(-1.2111726) q[3];
sx q[3];
rz(-2.6611879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2368161) q[0];
sx q[0];
rz(-0.95229709) q[0];
sx q[0];
rz(-3.0425518) q[0];
rz(-2.1095443) q[1];
sx q[1];
rz(-2.2910304) q[1];
sx q[1];
rz(-1.7031857) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28171646) q[0];
sx q[0];
rz(-1.8252354) q[0];
sx q[0];
rz(0.86812015) q[0];
rz(-pi) q[1];
rz(1.4969669) q[2];
sx q[2];
rz(-1.2880688) q[2];
sx q[2];
rz(1.7746995) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3308823) q[1];
sx q[1];
rz(-1.5473978) q[1];
sx q[1];
rz(0.7903654) q[1];
rz(-pi) q[2];
rz(-0.37495965) q[3];
sx q[3];
rz(-0.31878933) q[3];
sx q[3];
rz(-0.44170435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9859887) q[2];
sx q[2];
rz(-2.6376548) q[2];
sx q[2];
rz(-2.3206553) q[2];
rz(-1.4486897) q[3];
sx q[3];
rz(-2.4235348) q[3];
sx q[3];
rz(1.6528355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89707017) q[0];
sx q[0];
rz(-2.2269766) q[0];
sx q[0];
rz(-0.50265092) q[0];
rz(1.2333168) q[1];
sx q[1];
rz(-2.2063467) q[1];
sx q[1];
rz(0.9800235) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21979688) q[0];
sx q[0];
rz(-1.8077407) q[0];
sx q[0];
rz(-1.3150929) q[0];
rz(-0.29666846) q[2];
sx q[2];
rz(-0.7985332) q[2];
sx q[2];
rz(-2.1615596) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7240643) q[1];
sx q[1];
rz(-1.8019946) q[1];
sx q[1];
rz(2.3910012) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3408889) q[3];
sx q[3];
rz(-2.3977444) q[3];
sx q[3];
rz(-0.4901674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.22214733) q[2];
sx q[2];
rz(-1.28747) q[2];
sx q[2];
rz(1.7670828) q[2];
rz(1.8253271) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(-0.98178664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027243622) q[0];
sx q[0];
rz(-2.7482432) q[0];
sx q[0];
rz(2.223176) q[0];
rz(-0.20485993) q[1];
sx q[1];
rz(-2.0826976) q[1];
sx q[1];
rz(1.6580261) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433454) q[0];
sx q[0];
rz(-1.7172617) q[0];
sx q[0];
rz(-1.8634169) q[0];
x q[1];
rz(3.0181418) q[2];
sx q[2];
rz(-2.3268893) q[2];
sx q[2];
rz(-1.3208117) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10459722) q[1];
sx q[1];
rz(-0.87338398) q[1];
sx q[1];
rz(-3.0265977) q[1];
rz(-pi) q[2];
x q[2];
rz(0.056413944) q[3];
sx q[3];
rz(-2.4216363) q[3];
sx q[3];
rz(2.3145793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.001699) q[2];
sx q[2];
rz(-0.49109444) q[2];
sx q[2];
rz(1.3341382) q[2];
rz(-2.259518) q[3];
sx q[3];
rz(-0.69552723) q[3];
sx q[3];
rz(-1.2923366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081721574) q[0];
sx q[0];
rz(-2.7099755) q[0];
sx q[0];
rz(-0.01817848) q[0];
rz(-0.40953088) q[1];
sx q[1];
rz(-1.7117701) q[1];
sx q[1];
rz(0.10083625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026196711) q[0];
sx q[0];
rz(-1.5362722) q[0];
sx q[0];
rz(2.5071397) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42068215) q[2];
sx q[2];
rz(-1.5245588) q[2];
sx q[2];
rz(1.7939523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8162569) q[1];
sx q[1];
rz(-0.44202572) q[1];
sx q[1];
rz(-0.97953743) q[1];
x q[2];
rz(-2.0193897) q[3];
sx q[3];
rz(-0.27315419) q[3];
sx q[3];
rz(-2.7970683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8044901) q[2];
sx q[2];
rz(-0.10919315) q[2];
sx q[2];
rz(1.2479372) q[2];
rz(-1.2248056) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(0.065751806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.2357904) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(0.32050785) q[0];
rz(-0.39101741) q[1];
sx q[1];
rz(-1.1415488) q[1];
sx q[1];
rz(1.5919707) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3282916) q[0];
sx q[0];
rz(-1.3718954) q[0];
sx q[0];
rz(-1.3764253) q[0];
rz(-pi) q[1];
rz(0.28279998) q[2];
sx q[2];
rz(-0.53164266) q[2];
sx q[2];
rz(2.799019) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.446881) q[1];
sx q[1];
rz(-1.9384192) q[1];
sx q[1];
rz(3.047961) q[1];
rz(-pi) q[2];
rz(-3.1248766) q[3];
sx q[3];
rz(-2.7191945) q[3];
sx q[3];
rz(-1.0836386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.090342) q[2];
sx q[2];
rz(-1.0498468) q[2];
sx q[2];
rz(-2.7412097) q[2];
rz(2.9054437) q[3];
sx q[3];
rz(-3.0381687) q[3];
sx q[3];
rz(-1.0108488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.878933) q[0];
sx q[0];
rz(-0.50983179) q[0];
sx q[0];
rz(-1.8806993) q[0];
rz(-0.48514584) q[1];
sx q[1];
rz(-0.79882516) q[1];
sx q[1];
rz(-0.47732236) q[1];
rz(3.015949) q[2];
sx q[2];
rz(-1.7775272) q[2];
sx q[2];
rz(0.20128332) q[2];
rz(0.063592576) q[3];
sx q[3];
rz(-1.7983754) q[3];
sx q[3];
rz(1.6354431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

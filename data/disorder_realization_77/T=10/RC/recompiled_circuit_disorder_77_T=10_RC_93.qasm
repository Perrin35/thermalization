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
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6457155) q[0];
sx q[0];
rz(-2.0464532) q[0];
sx q[0];
rz(0.23947421) q[0];
rz(-2.9287455) q[2];
sx q[2];
rz(-0.93570645) q[2];
sx q[2];
rz(-1.1320621) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1311156) q[1];
sx q[1];
rz(-1.3033723) q[1];
sx q[1];
rz(-2.1285776) q[1];
rz(-1.6151186) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(-2.0093902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7636259) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(0.4326694) q[2];
rz(1.1928605) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6137961) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(-1.312785) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.1516494) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0002999) q[0];
sx q[0];
rz(-2.2465758) q[0];
sx q[0];
rz(-2.7594901) q[0];
rz(-pi) q[1];
rz(-0.58354124) q[2];
sx q[2];
rz(-1.984664) q[2];
sx q[2];
rz(-1.5838503) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18560219) q[1];
sx q[1];
rz(-2.5163109) q[1];
sx q[1];
rz(0.76693265) q[1];
rz(-pi) q[2];
rz(-1.1851951) q[3];
sx q[3];
rz(-1.10023) q[3];
sx q[3];
rz(-0.040599559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1318704) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-0.22182626) q[2];
rz(0.37718537) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(-2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(-0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(3.085014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6037613) q[0];
sx q[0];
rz(-1.0147525) q[0];
sx q[0];
rz(-2.6105196) q[0];
rz(-pi) q[1];
rz(-1.7043731) q[2];
sx q[2];
rz(-1.3165511) q[2];
sx q[2];
rz(0.1453407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6371582) q[1];
sx q[1];
rz(-1.27099) q[1];
sx q[1];
rz(-1.7791041) q[1];
rz(-1.9583086) q[3];
sx q[3];
rz(-2.1790677) q[3];
sx q[3];
rz(-1.1063948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.27292192) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(0.92612129) q[2];
rz(0.55666322) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(-2.0986957) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(0.74209374) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-0.46359584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.050042) q[0];
sx q[0];
rz(-2.3766962) q[0];
sx q[0];
rz(1.1372304) q[0];
x q[1];
rz(-1.5092588) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(-2.5663944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3568748) q[1];
sx q[1];
rz(-0.80662913) q[1];
sx q[1];
rz(-2.9033317) q[1];
x q[2];
rz(0.87426825) q[3];
sx q[3];
rz(-0.6797176) q[3];
sx q[3];
rz(-3.0391039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7745557) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(-3.0130623) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(-0.29770011) q[0];
rz(-0.4822576) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(0.94435) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1221065) q[0];
sx q[0];
rz(-1.6158316) q[0];
sx q[0];
rz(-1.3615863) q[0];
x q[1];
rz(-1.1733426) q[2];
sx q[2];
rz(-0.83565088) q[2];
sx q[2];
rz(-0.022692516) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.797232) q[1];
sx q[1];
rz(-1.6318983) q[1];
sx q[1];
rz(-3.1394715) q[1];
x q[2];
rz(-0.43290187) q[3];
sx q[3];
rz(-1.8042943) q[3];
sx q[3];
rz(1.3732861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8828316) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(3.1392858) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(0.58445245) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-3.086673) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7391101) q[0];
sx q[0];
rz(-1.6533028) q[0];
sx q[0];
rz(1.301618) q[0];
rz(-pi) q[1];
rz(1.8940077) q[2];
sx q[2];
rz(-1.5885457) q[2];
sx q[2];
rz(0.36349597) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1850486) q[1];
sx q[1];
rz(-1.1951606) q[1];
sx q[1];
rz(-0.59021414) q[1];
x q[2];
rz(2.1543598) q[3];
sx q[3];
rz(-1.9658486) q[3];
sx q[3];
rz(0.46056718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75446689) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(-2.1248655) q[2];
rz(0.54404849) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465437) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(0.28453919) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(2.231266) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9757662) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(1.5058917) q[0];
rz(2.8001627) q[2];
sx q[2];
rz(-1.9651946) q[2];
sx q[2];
rz(-0.1614801) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6701723) q[1];
sx q[1];
rz(-2.6176665) q[1];
sx q[1];
rz(2.7736204) q[1];
rz(-pi) q[2];
rz(2.4207553) q[3];
sx q[3];
rz(-1.1158873) q[3];
sx q[3];
rz(2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(2.2195623) q[2];
rz(-0.56728029) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(-1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(-2.5091876) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(0.30050373) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104615) q[0];
sx q[0];
rz(-2.3346402) q[0];
sx q[0];
rz(2.0956844) q[0];
x q[1];
rz(-2.3807316) q[2];
sx q[2];
rz(-1.0957452) q[2];
sx q[2];
rz(2.8617815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.955519) q[1];
sx q[1];
rz(-1.7454073) q[1];
sx q[1];
rz(-0.34592918) q[1];
rz(-pi) q[2];
rz(0.31605966) q[3];
sx q[3];
rz(-0.83871597) q[3];
sx q[3];
rz(2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5552716) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-2.3596181) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5732116) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-2.382747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0046783) q[0];
sx q[0];
rz(-1.5877962) q[0];
sx q[0];
rz(-1.9949811) q[0];
x q[1];
rz(-0.27300948) q[2];
sx q[2];
rz(-2.5148279) q[2];
sx q[2];
rz(3.0221456) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85941852) q[1];
sx q[1];
rz(-2.4255883) q[1];
sx q[1];
rz(2.9031258) q[1];
rz(1.5893448) q[3];
sx q[3];
rz(-2.4759001) q[3];
sx q[3];
rz(-1.9513643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(2.0786044) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(0.075335659) q[0];
rz(0.8967337) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(-2.5316701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2595554) q[0];
sx q[0];
rz(-1.8543188) q[0];
sx q[0];
rz(-0.77772227) q[0];
x q[1];
rz(-0.62060771) q[2];
sx q[2];
rz(-2.6891516) q[2];
sx q[2];
rz(2.0394182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1249485) q[1];
sx q[1];
rz(-2.2955107) q[1];
sx q[1];
rz(-1.8925341) q[1];
rz(-pi) q[2];
x q[2];
rz(0.742357) q[3];
sx q[3];
rz(-1.1872113) q[3];
sx q[3];
rz(2.50768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23218368) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(-2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80355766) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(2.4907885) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(1.4177633) q[2];
sx q[2];
rz(-2.9433708) q[2];
sx q[2];
rz(2.3797258) q[2];
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

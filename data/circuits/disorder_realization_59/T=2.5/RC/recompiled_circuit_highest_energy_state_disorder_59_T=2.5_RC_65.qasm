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
rz(2.4177999) q[0];
sx q[0];
rz(-1.1136709) q[0];
sx q[0];
rz(-0.96473515) q[0];
rz(-6.3568883) q[1];
sx q[1];
rz(2.5001723) q[1];
sx q[1];
rz(14.060798) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.651123) q[0];
sx q[0];
rz(-2.2283366) q[0];
sx q[0];
rz(1.0762908) q[0];
x q[1];
rz(-1.2624717) q[2];
sx q[2];
rz(-2.8092896) q[2];
sx q[2];
rz(0.52969709) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9656417) q[1];
sx q[1];
rz(-2.0276105) q[1];
sx q[1];
rz(-2.5547002) q[1];
rz(-pi) q[2];
x q[2];
rz(1.649664) q[3];
sx q[3];
rz(-0.53073629) q[3];
sx q[3];
rz(-2.4473878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4045495) q[2];
sx q[2];
rz(-1.5207542) q[2];
sx q[2];
rz(2.9250195) q[2];
rz(-0.51566044) q[3];
sx q[3];
rz(-2.0629081) q[3];
sx q[3];
rz(0.083958538) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11230042) q[0];
sx q[0];
rz(-0.071542112) q[0];
sx q[0];
rz(1.4805502) q[0];
rz(0.59313613) q[1];
sx q[1];
rz(-0.99855223) q[1];
sx q[1];
rz(0.28233972) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823989) q[0];
sx q[0];
rz(-1.6983573) q[0];
sx q[0];
rz(1.5414627) q[0];
rz(-pi) q[1];
rz(-0.13377269) q[2];
sx q[2];
rz(-1.7672799) q[2];
sx q[2];
rz(-1.5143732) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8806539) q[1];
sx q[1];
rz(-1.0323413) q[1];
sx q[1];
rz(1.4101342) q[1];
rz(-pi) q[2];
rz(2.2766791) q[3];
sx q[3];
rz(-2.5967715) q[3];
sx q[3];
rz(-0.13710216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7181519) q[2];
sx q[2];
rz(-0.56090063) q[2];
sx q[2];
rz(-0.94914985) q[2];
rz(-0.078350457) q[3];
sx q[3];
rz(-1.4698942) q[3];
sx q[3];
rz(1.2300434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3345728) q[0];
sx q[0];
rz(-0.0581352) q[0];
sx q[0];
rz(0.19244254) q[0];
rz(2.1678534) q[1];
sx q[1];
rz(-2.1111919) q[1];
sx q[1];
rz(1.8691501) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26548037) q[0];
sx q[0];
rz(-2.2372549) q[0];
sx q[0];
rz(-1.8811532) q[0];
rz(-pi) q[1];
rz(-0.22120938) q[2];
sx q[2];
rz(-2.2941046) q[2];
sx q[2];
rz(2.4257367) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7682926) q[1];
sx q[1];
rz(-1.321861) q[1];
sx q[1];
rz(1.8087923) q[1];
rz(3.104464) q[3];
sx q[3];
rz(-1.1151966) q[3];
sx q[3];
rz(1.451481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56866139) q[2];
sx q[2];
rz(-2.2243786) q[2];
sx q[2];
rz(0.60079637) q[2];
rz(0.83493799) q[3];
sx q[3];
rz(-0.91885126) q[3];
sx q[3];
rz(0.44343534) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610953) q[0];
sx q[0];
rz(-2.0090071) q[0];
sx q[0];
rz(1.7219211) q[0];
rz(2.8555866) q[1];
sx q[1];
rz(-1.1347457) q[1];
sx q[1];
rz(2.6223415) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70066455) q[0];
sx q[0];
rz(-1.6608149) q[0];
sx q[0];
rz(0.16914455) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18188484) q[2];
sx q[2];
rz(-1.4684235) q[2];
sx q[2];
rz(-0.052396862) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45471474) q[1];
sx q[1];
rz(-1.5022394) q[1];
sx q[1];
rz(2.6297533) q[1];
rz(-pi) q[2];
x q[2];
rz(0.052922225) q[3];
sx q[3];
rz(-0.54659802) q[3];
sx q[3];
rz(-1.6569759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9754192) q[2];
sx q[2];
rz(-0.68048802) q[2];
sx q[2];
rz(0.55230459) q[2];
rz(2.35899) q[3];
sx q[3];
rz(-1.3026594) q[3];
sx q[3];
rz(2.9569614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85823378) q[0];
sx q[0];
rz(-2.5193546) q[0];
sx q[0];
rz(1.8587814) q[0];
rz(-1.2239617) q[1];
sx q[1];
rz(-1.6061648) q[1];
sx q[1];
rz(0.70972365) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3027355) q[0];
sx q[0];
rz(-1.5488911) q[0];
sx q[0];
rz(3.1126853) q[0];
rz(-2.6733396) q[2];
sx q[2];
rz(-1.5522458) q[2];
sx q[2];
rz(1.5520688) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9319871) q[1];
sx q[1];
rz(-2.8199204) q[1];
sx q[1];
rz(-1.7965921) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.157326) q[3];
sx q[3];
rz(-1.4063104) q[3];
sx q[3];
rz(2.2604347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84511406) q[2];
sx q[2];
rz(-0.60135403) q[2];
sx q[2];
rz(2.5884886) q[2];
rz(0.84135711) q[3];
sx q[3];
rz(-1.0480806) q[3];
sx q[3];
rz(1.1355737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3784921) q[0];
sx q[0];
rz(-1.0569514) q[0];
sx q[0];
rz(-2.3665449) q[0];
rz(1.4729602) q[1];
sx q[1];
rz(-2.5170363) q[1];
sx q[1];
rz(-0.83736173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47507163) q[0];
sx q[0];
rz(-1.6672264) q[0];
sx q[0];
rz(0.94397105) q[0];
rz(2.7816394) q[2];
sx q[2];
rz(-2.6591566) q[2];
sx q[2];
rz(-3.0878518) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7720799) q[1];
sx q[1];
rz(-1.9926148) q[1];
sx q[1];
rz(-1.1320516) q[1];
rz(1.9894137) q[3];
sx q[3];
rz(-1.3411203) q[3];
sx q[3];
rz(2.9269689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.9019292) q[2];
sx q[2];
rz(-1.1666802) q[2];
sx q[2];
rz(-0.76266328) q[2];
rz(2.934382) q[3];
sx q[3];
rz(-1.8800294) q[3];
sx q[3];
rz(0.60780779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36829456) q[0];
sx q[0];
rz(-0.93450707) q[0];
sx q[0];
rz(1.4542798) q[0];
rz(1.2794718) q[1];
sx q[1];
rz(-2.2984633) q[1];
sx q[1];
rz(-1.4131193) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18943809) q[0];
sx q[0];
rz(-1.3690152) q[0];
sx q[0];
rz(0.39315572) q[0];
rz(-pi) q[1];
rz(-2.3033875) q[2];
sx q[2];
rz(-1.1229441) q[2];
sx q[2];
rz(-3.107576) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5190025) q[1];
sx q[1];
rz(-1.1347924) q[1];
sx q[1];
rz(1.9080129) q[1];
rz(-2.0094447) q[3];
sx q[3];
rz(-1.2474212) q[3];
sx q[3];
rz(-1.974859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9965245) q[2];
sx q[2];
rz(-0.86864305) q[2];
sx q[2];
rz(0.16150148) q[2];
rz(-0.7343556) q[3];
sx q[3];
rz(-2.2752094) q[3];
sx q[3];
rz(1.3041147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2187918) q[0];
sx q[0];
rz(-0.53023338) q[0];
sx q[0];
rz(2.9042322) q[0];
rz(-0.78422076) q[1];
sx q[1];
rz(-1.3619245) q[1];
sx q[1];
rz(1.721419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7331384) q[0];
sx q[0];
rz(-1.2722005) q[0];
sx q[0];
rz(-0.51513735) q[0];
rz(2.2591822) q[2];
sx q[2];
rz(-1.3974481) q[2];
sx q[2];
rz(-1.8515996) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.776122) q[1];
sx q[1];
rz(-2.104722) q[1];
sx q[1];
rz(-1.9129685) q[1];
x q[2];
rz(-2.0886252) q[3];
sx q[3];
rz(-1.1204567) q[3];
sx q[3];
rz(-0.88847775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4935442) q[2];
sx q[2];
rz(-1.8914787) q[2];
sx q[2];
rz(1.9165215) q[2];
rz(-1.8152292) q[3];
sx q[3];
rz(-2.1882961) q[3];
sx q[3];
rz(1.1939322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61106435) q[0];
sx q[0];
rz(-3.0623797) q[0];
sx q[0];
rz(0.63842574) q[0];
rz(3.0910659) q[1];
sx q[1];
rz(-1.9332998) q[1];
sx q[1];
rz(2.5343177) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3359708) q[0];
sx q[0];
rz(-2.6329649) q[0];
sx q[0];
rz(-0.87368272) q[0];
rz(-pi) q[1];
rz(0.88596099) q[2];
sx q[2];
rz(-1.6185624) q[2];
sx q[2];
rz(1.1807962) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7672459) q[1];
sx q[1];
rz(-2.8082529) q[1];
sx q[1];
rz(-1.7872168) q[1];
x q[2];
rz(2.3361198) q[3];
sx q[3];
rz(-2.8777524) q[3];
sx q[3];
rz(-1.4782617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4813469) q[2];
sx q[2];
rz(-1.1911743) q[2];
sx q[2];
rz(-2.353239) q[2];
rz(-2.3039019) q[3];
sx q[3];
rz(-2.3205784) q[3];
sx q[3];
rz(-1.2825509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4072708) q[0];
sx q[0];
rz(-2.4096074) q[0];
sx q[0];
rz(0.5087854) q[0];
rz(1.5599627) q[1];
sx q[1];
rz(-2.1211233) q[1];
sx q[1];
rz(2.7324953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74241766) q[0];
sx q[0];
rz(-1.3361201) q[0];
sx q[0];
rz(0.15813903) q[0];
x q[1];
rz(0.47637911) q[2];
sx q[2];
rz(-2.0198108) q[2];
sx q[2];
rz(0.19461122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.075114014) q[1];
sx q[1];
rz(-1.7914322) q[1];
sx q[1];
rz(-1.5815939) q[1];
rz(-pi) q[2];
rz(-0.76830058) q[3];
sx q[3];
rz(-2.6242206) q[3];
sx q[3];
rz(-0.20341541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.540655) q[2];
sx q[2];
rz(-1.6161852) q[2];
sx q[2];
rz(1.8589004) q[2];
rz(-0.013785275) q[3];
sx q[3];
rz(-1.1252334) q[3];
sx q[3];
rz(1.8085326) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8448821) q[0];
sx q[0];
rz(-1.7322576) q[0];
sx q[0];
rz(-2.5175293) q[0];
rz(0.86112549) q[1];
sx q[1];
rz(-0.54665165) q[1];
sx q[1];
rz(0.13269592) q[1];
rz(-1.4693946) q[2];
sx q[2];
rz(-1.8894926) q[2];
sx q[2];
rz(2.8429902) q[2];
rz(0.63152159) q[3];
sx q[3];
rz(-2.1248264) q[3];
sx q[3];
rz(-0.2223224) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

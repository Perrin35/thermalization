OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3026128) q[0];
sx q[0];
rz(-1.4361199) q[0];
sx q[0];
rz(0.21685313) q[0];
rz(-0.4878374) q[1];
sx q[1];
rz(-0.87895972) q[1];
sx q[1];
rz(2.9882957) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3756806) q[0];
sx q[0];
rz(-1.8366251) q[0];
sx q[0];
rz(1.2714809) q[0];
rz(0.4699444) q[2];
sx q[2];
rz(-2.0389839) q[2];
sx q[2];
rz(2.9681457) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.98521122) q[1];
sx q[1];
rz(-2.0655246) q[1];
sx q[1];
rz(-1.8534766) q[1];
rz(-0.1568832) q[3];
sx q[3];
rz(-2.3290754) q[3];
sx q[3];
rz(2.7650327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45120254) q[2];
sx q[2];
rz(-2.5610552) q[2];
sx q[2];
rz(2.2817877) q[2];
rz(2.9016923) q[3];
sx q[3];
rz(-0.71687117) q[3];
sx q[3];
rz(-0.2828323) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4848223) q[0];
sx q[0];
rz(-1.7253933) q[0];
sx q[0];
rz(-2.2221478) q[0];
rz(0.8473618) q[1];
sx q[1];
rz(-0.58440009) q[1];
sx q[1];
rz(-1.1712317) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4019415) q[0];
sx q[0];
rz(-1.8293132) q[0];
sx q[0];
rz(3.0895732) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33185766) q[2];
sx q[2];
rz(-2.2584791) q[2];
sx q[2];
rz(0.0092384641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6379063) q[1];
sx q[1];
rz(-1.827938) q[1];
sx q[1];
rz(1.2127582) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1055206) q[3];
sx q[3];
rz(-2.6031446) q[3];
sx q[3];
rz(-3.127436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.072711572) q[2];
sx q[2];
rz(-0.88901192) q[2];
sx q[2];
rz(1.2129126) q[2];
rz(-2.9291901) q[3];
sx q[3];
rz(-2.0272777) q[3];
sx q[3];
rz(-0.67306486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251959) q[0];
sx q[0];
rz(-1.9028417) q[0];
sx q[0];
rz(-2.5585001) q[0];
rz(-1.8816226) q[1];
sx q[1];
rz(-0.59695736) q[1];
sx q[1];
rz(-1.3139542) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010600032) q[0];
sx q[0];
rz(-0.48904648) q[0];
sx q[0];
rz(-0.84757342) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.407786) q[2];
sx q[2];
rz(-1.5948405) q[2];
sx q[2];
rz(0.95685416) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3594234) q[1];
sx q[1];
rz(-0.12453989) q[1];
sx q[1];
rz(1.4119536) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2168808) q[3];
sx q[3];
rz(-2.6887213) q[3];
sx q[3];
rz(2.4726506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6446357) q[2];
sx q[2];
rz(-2.3778215) q[2];
sx q[2];
rz(2.606707) q[2];
rz(-2.6584117) q[3];
sx q[3];
rz(-1.6202241) q[3];
sx q[3];
rz(0.11972891) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0508761) q[0];
sx q[0];
rz(-3.0829939) q[0];
sx q[0];
rz(1.6706049) q[0];
rz(-1.927467) q[1];
sx q[1];
rz(-1.7045226) q[1];
sx q[1];
rz(2.6953183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4497919) q[0];
sx q[0];
rz(-1.746382) q[0];
sx q[0];
rz(3.0719766) q[0];
rz(0.52881119) q[2];
sx q[2];
rz(-1.5545648) q[2];
sx q[2];
rz(-1.123614) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.41149263) q[1];
sx q[1];
rz(-0.92918438) q[1];
sx q[1];
rz(1.8909251) q[1];
rz(1.5851444) q[3];
sx q[3];
rz(-1.8135462) q[3];
sx q[3];
rz(-1.2402463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92397592) q[2];
sx q[2];
rz(-2.6805704) q[2];
sx q[2];
rz(2.8154742) q[2];
rz(2.7255132) q[3];
sx q[3];
rz(-1.5030428) q[3];
sx q[3];
rz(-2.3638341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.07051) q[0];
sx q[0];
rz(-1.5639045) q[0];
sx q[0];
rz(-0.34410205) q[0];
rz(0.35585078) q[1];
sx q[1];
rz(-0.75332037) q[1];
sx q[1];
rz(-2.8867302) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0094285) q[0];
sx q[0];
rz(-0.71064204) q[0];
sx q[0];
rz(-0.406213) q[0];
rz(-1.1084313) q[2];
sx q[2];
rz(-1.6605262) q[2];
sx q[2];
rz(-1.7347908) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5755467) q[1];
sx q[1];
rz(-2.0064163) q[1];
sx q[1];
rz(0.91710414) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1625103) q[3];
sx q[3];
rz(-2.2934543) q[3];
sx q[3];
rz(1.6734139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7039589) q[2];
sx q[2];
rz(-0.86486977) q[2];
sx q[2];
rz(3.051009) q[2];
rz(-0.80638742) q[3];
sx q[3];
rz(-1.3017637) q[3];
sx q[3];
rz(1.8346132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059747132) q[0];
sx q[0];
rz(-1.9192001) q[0];
sx q[0];
rz(2.5229689) q[0];
rz(0.6126569) q[1];
sx q[1];
rz(-2.5109992) q[1];
sx q[1];
rz(2.9887066) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3090749) q[0];
sx q[0];
rz(-0.33402696) q[0];
sx q[0];
rz(-0.79132737) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15494187) q[2];
sx q[2];
rz(-2.5563588) q[2];
sx q[2];
rz(-2.3335148) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5093435) q[1];
sx q[1];
rz(-1.9652838) q[1];
sx q[1];
rz(-0.057827397) q[1];
rz(2.1508407) q[3];
sx q[3];
rz(-1.5501805) q[3];
sx q[3];
rz(-1.0643268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20554081) q[2];
sx q[2];
rz(-2.3369393) q[2];
sx q[2];
rz(-1.12977) q[2];
rz(1.0063082) q[3];
sx q[3];
rz(-1.8215424) q[3];
sx q[3];
rz(1.6442851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048700843) q[0];
sx q[0];
rz(-2.0755656) q[0];
sx q[0];
rz(-0.14376465) q[0];
rz(0.20214209) q[1];
sx q[1];
rz(-0.527924) q[1];
sx q[1];
rz(-1.082083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7162299) q[0];
sx q[0];
rz(-0.83304616) q[0];
sx q[0];
rz(-2.0701755) q[0];
rz(0.80777709) q[2];
sx q[2];
rz(-0.91726724) q[2];
sx q[2];
rz(-1.8612922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.52817594) q[1];
sx q[1];
rz(-0.74494637) q[1];
sx q[1];
rz(1.2553535) q[1];
x q[2];
rz(2.4025926) q[3];
sx q[3];
rz(-2.0428052) q[3];
sx q[3];
rz(-1.3907741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9988592) q[2];
sx q[2];
rz(-1.1799246) q[2];
sx q[2];
rz(-0.72706968) q[2];
rz(-1.8266228) q[3];
sx q[3];
rz(-1.246779) q[3];
sx q[3];
rz(-1.6453913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24406782) q[0];
sx q[0];
rz(-2.0638564) q[0];
sx q[0];
rz(-1.1429853) q[0];
rz(2.9179528) q[1];
sx q[1];
rz(-0.63839212) q[1];
sx q[1];
rz(-1.2248096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6055034) q[0];
sx q[0];
rz(-1.8538626) q[0];
sx q[0];
rz(-0.27293954) q[0];
rz(-1.9838422) q[2];
sx q[2];
rz(-0.89282477) q[2];
sx q[2];
rz(1.7630446) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.599419) q[1];
sx q[1];
rz(-1.6288593) q[1];
sx q[1];
rz(1.7752613) q[1];
rz(-pi) q[2];
rz(0.16815987) q[3];
sx q[3];
rz(-2.5736397) q[3];
sx q[3];
rz(2.0310543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58330047) q[2];
sx q[2];
rz(-1.4536828) q[2];
sx q[2];
rz(0.80201403) q[2];
rz(-1.924104) q[3];
sx q[3];
rz(-0.13974443) q[3];
sx q[3];
rz(-1.2453992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8121346) q[0];
sx q[0];
rz(-1.7725002) q[0];
sx q[0];
rz(1.580397) q[0];
rz(2.2691057) q[1];
sx q[1];
rz(-0.28631887) q[1];
sx q[1];
rz(-2.7365541) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0035343) q[0];
sx q[0];
rz(-2.0857607) q[0];
sx q[0];
rz(-2.1206417) q[0];
rz(-2.9895859) q[2];
sx q[2];
rz(-0.59796732) q[2];
sx q[2];
rz(2.1961308) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.424347) q[1];
sx q[1];
rz(-0.97386347) q[1];
sx q[1];
rz(-2.2655377) q[1];
x q[2];
rz(-1.551335) q[3];
sx q[3];
rz(-1.4238589) q[3];
sx q[3];
rz(-2.9732239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.27434719) q[2];
sx q[2];
rz(-0.87928191) q[2];
sx q[2];
rz(0.61152968) q[2];
rz(-1.3036171) q[3];
sx q[3];
rz(-2.7795064) q[3];
sx q[3];
rz(-1.1360137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63854727) q[0];
sx q[0];
rz(-1.2940116) q[0];
sx q[0];
rz(3.0882623) q[0];
rz(-2.5697925) q[1];
sx q[1];
rz(-1.6245533) q[1];
sx q[1];
rz(1.8624051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17634468) q[0];
sx q[0];
rz(-1.2851479) q[0];
sx q[0];
rz(-1.833451) q[0];
rz(-pi) q[1];
rz(-1.2917742) q[2];
sx q[2];
rz(-2.6447562) q[2];
sx q[2];
rz(1.9328839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0030583) q[1];
sx q[1];
rz(-1.6489677) q[1];
sx q[1];
rz(-0.98346904) q[1];
rz(-2.0883858) q[3];
sx q[3];
rz(-2.7206051) q[3];
sx q[3];
rz(1.8618575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6751487) q[2];
sx q[2];
rz(-1.6089336) q[2];
sx q[2];
rz(-0.68501985) q[2];
rz(0.78392616) q[3];
sx q[3];
rz(-2.7214366) q[3];
sx q[3];
rz(2.4043731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53981275) q[0];
sx q[0];
rz(-1.89986) q[0];
sx q[0];
rz(1.4358406) q[0];
rz(2.336179) q[1];
sx q[1];
rz(-0.46331159) q[1];
sx q[1];
rz(0.70809271) q[1];
rz(1.3546181) q[2];
sx q[2];
rz(-1.8541649) q[2];
sx q[2];
rz(-0.46021663) q[2];
rz(0.2652771) q[3];
sx q[3];
rz(-0.42828413) q[3];
sx q[3];
rz(-0.47676906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

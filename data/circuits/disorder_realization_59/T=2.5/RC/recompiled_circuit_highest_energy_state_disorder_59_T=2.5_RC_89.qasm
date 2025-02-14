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
rz(2.1768575) q[0];
rz(-6.3568883) q[1];
sx q[1];
rz(2.5001723) q[1];
sx q[1];
rz(14.060798) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.651123) q[0];
sx q[0];
rz(-0.91325608) q[0];
sx q[0];
rz(-1.0762908) q[0];
rz(-pi) q[1];
x q[1];
rz(3.037248) q[2];
sx q[2];
rz(-1.2547224) q[2];
sx q[2];
rz(2.2869273) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19077483) q[1];
sx q[1];
rz(-0.72682805) q[1];
sx q[1];
rz(0.72587691) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0414181) q[3];
sx q[3];
rz(-1.6106859) q[3];
sx q[3];
rz(-2.1969469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4045495) q[2];
sx q[2];
rz(-1.5207542) q[2];
sx q[2];
rz(-0.21657319) q[2];
rz(0.51566044) q[3];
sx q[3];
rz(-1.0786846) q[3];
sx q[3];
rz(0.083958538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0292922) q[0];
sx q[0];
rz(-3.0700505) q[0];
sx q[0];
rz(-1.4805502) q[0];
rz(-2.5484565) q[1];
sx q[1];
rz(-2.1430404) q[1];
sx q[1];
rz(-0.28233972) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7337228) q[0];
sx q[0];
rz(-1.5417011) q[0];
sx q[0];
rz(-0.12761527) q[0];
rz(-1.3725874) q[2];
sx q[2];
rz(-1.4396126) q[2];
sx q[2];
rz(3.1114357) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9548492) q[1];
sx q[1];
rz(-2.5819518) q[1];
sx q[1];
rz(-0.26167027) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1386747) q[3];
sx q[3];
rz(-1.2279155) q[3];
sx q[3];
rz(-2.0636095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4234408) q[2];
sx q[2];
rz(-2.580692) q[2];
sx q[2];
rz(-0.94914985) q[2];
rz(-3.0632422) q[3];
sx q[3];
rz(-1.6716985) q[3];
sx q[3];
rz(-1.9115492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80701989) q[0];
sx q[0];
rz(-0.0581352) q[0];
sx q[0];
rz(-2.9491501) q[0];
rz(-2.1678534) q[1];
sx q[1];
rz(-1.0304008) q[1];
sx q[1];
rz(-1.2724426) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26548037) q[0];
sx q[0];
rz(-0.9043378) q[0];
sx q[0];
rz(-1.2604395) q[0];
rz(-0.22120938) q[2];
sx q[2];
rz(-2.2941046) q[2];
sx q[2];
rz(2.4257367) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0037937) q[1];
sx q[1];
rz(-1.8013211) q[1];
sx q[1];
rz(-2.8857438) q[1];
x q[2];
rz(-0.037128673) q[3];
sx q[3];
rz(-2.026396) q[3];
sx q[3];
rz(1.6901117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5729313) q[2];
sx q[2];
rz(-0.91721407) q[2];
sx q[2];
rz(2.5407963) q[2];
rz(2.3066547) q[3];
sx q[3];
rz(-2.2227414) q[3];
sx q[3];
rz(-2.6981573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610953) q[0];
sx q[0];
rz(-2.0090071) q[0];
sx q[0];
rz(-1.4196716) q[0];
rz(-2.8555866) q[1];
sx q[1];
rz(-2.0068469) q[1];
sx q[1];
rz(2.6223415) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38558593) q[0];
sx q[0];
rz(-2.950188) q[0];
sx q[0];
rz(-2.6494104) q[0];
rz(-0.18188484) q[2];
sx q[2];
rz(-1.6731691) q[2];
sx q[2];
rz(-3.0891958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6868779) q[1];
sx q[1];
rz(-1.5022394) q[1];
sx q[1];
rz(2.6297533) q[1];
rz(-1.5386228) q[3];
sx q[3];
rz(-2.1165427) q[3];
sx q[3];
rz(1.4226899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.1661735) q[2];
sx q[2];
rz(-0.68048802) q[2];
sx q[2];
rz(-2.5892881) q[2];
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
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85823378) q[0];
sx q[0];
rz(-0.62223804) q[0];
sx q[0];
rz(1.8587814) q[0];
rz(1.917631) q[1];
sx q[1];
rz(-1.5354278) q[1];
sx q[1];
rz(-0.70972365) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8741651) q[0];
sx q[0];
rz(-1.5418959) q[0];
sx q[0];
rz(-1.5927107) q[0];
x q[1];
rz(-0.041083606) q[2];
sx q[2];
rz(-2.6729995) q[2];
sx q[2];
rz(0.017931708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2096055) q[1];
sx q[1];
rz(-2.8199204) q[1];
sx q[1];
rz(-1.7965921) q[1];
rz(-pi) q[2];
rz(0.17931314) q[3];
sx q[3];
rz(-1.1632405) q[3];
sx q[3];
rz(2.5236764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84511406) q[2];
sx q[2];
rz(-2.5402386) q[2];
sx q[2];
rz(0.55310407) q[2];
rz(-0.84135711) q[3];
sx q[3];
rz(-1.0480806) q[3];
sx q[3];
rz(2.0060189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(2.3784921) q[0];
sx q[0];
rz(-1.0569514) q[0];
sx q[0];
rz(-2.3665449) q[0];
rz(-1.4729602) q[1];
sx q[1];
rz(-2.5170363) q[1];
sx q[1];
rz(0.83736173) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178029) q[0];
sx q[0];
rz(-0.63321165) q[0];
sx q[0];
rz(1.4073611) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35995324) q[2];
sx q[2];
rz(-0.48243603) q[2];
sx q[2];
rz(-3.0878518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0114767) q[1];
sx q[1];
rz(-1.9687879) q[1];
sx q[1];
rz(2.6813862) q[1];
x q[2];
rz(1.0488408) q[3];
sx q[3];
rz(-2.6673954) q[3];
sx q[3];
rz(1.8291345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2396635) q[2];
sx q[2];
rz(-1.9749125) q[2];
sx q[2];
rz(2.3789294) q[2];
rz(-0.20721063) q[3];
sx q[3];
rz(-1.8800294) q[3];
sx q[3];
rz(0.60780779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7732981) q[0];
sx q[0];
rz(-0.93450707) q[0];
sx q[0];
rz(1.6873129) q[0];
rz(-1.2794718) q[1];
sx q[1];
rz(-0.8431294) q[1];
sx q[1];
rz(1.7284733) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4642893) q[0];
sx q[0];
rz(-1.9555518) q[0];
sx q[0];
rz(1.3528525) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57369487) q[2];
sx q[2];
rz(-0.92373935) q[2];
sx q[2];
rz(1.1653125) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3403875) q[1];
sx q[1];
rz(-1.875393) q[1];
sx q[1];
rz(-2.6829863) q[1];
rz(2.0094447) q[3];
sx q[3];
rz(-1.2474212) q[3];
sx q[3];
rz(-1.1667337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9965245) q[2];
sx q[2];
rz(-2.2729496) q[2];
sx q[2];
rz(-0.16150148) q[2];
rz(0.7343556) q[3];
sx q[3];
rz(-2.2752094) q[3];
sx q[3];
rz(-1.3041147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9228009) q[0];
sx q[0];
rz(-2.6113593) q[0];
sx q[0];
rz(-2.9042322) q[0];
rz(-2.3573719) q[1];
sx q[1];
rz(-1.7796681) q[1];
sx q[1];
rz(-1.4201737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6415757) q[0];
sx q[0];
rz(-0.58860129) q[0];
sx q[0];
rz(0.55843784) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9186219) q[2];
sx q[2];
rz(-0.89467773) q[2];
sx q[2];
rz(-2.7198426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3846459) q[1];
sx q[1];
rz(-1.8638041) q[1];
sx q[1];
rz(-2.5811367) q[1];
rz(-2.3444207) q[3];
sx q[3];
rz(-2.4690921) q[3];
sx q[3];
rz(-0.029880015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64804849) q[2];
sx q[2];
rz(-1.8914787) q[2];
sx q[2];
rz(-1.9165215) q[2];
rz(-1.3263634) q[3];
sx q[3];
rz(-0.95329657) q[3];
sx q[3];
rz(-1.9476604) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61106435) q[0];
sx q[0];
rz(-3.0623797) q[0];
sx q[0];
rz(2.5031669) q[0];
rz(-3.0910659) q[1];
sx q[1];
rz(-1.9332998) q[1];
sx q[1];
rz(0.60727492) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3359708) q[0];
sx q[0];
rz(-2.6329649) q[0];
sx q[0];
rz(2.2679099) q[0];
x q[1];
rz(-1.6462244) q[2];
sx q[2];
rz(-0.68623073) q[2];
sx q[2];
rz(-2.8099912) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0083976204) q[1];
sx q[1];
rz(-1.6411157) q[1];
sx q[1];
rz(-1.2446872) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3361198) q[3];
sx q[3];
rz(-2.8777524) q[3];
sx q[3];
rz(1.4782617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4813469) q[2];
sx q[2];
rz(-1.1911743) q[2];
sx q[2];
rz(0.78835362) q[2];
rz(-0.83769074) q[3];
sx q[3];
rz(-0.82101429) q[3];
sx q[3];
rz(-1.2825509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4072708) q[0];
sx q[0];
rz(-0.73198524) q[0];
sx q[0];
rz(-0.5087854) q[0];
rz(-1.5599627) q[1];
sx q[1];
rz(-1.0204693) q[1];
sx q[1];
rz(2.7324953) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0003153) q[0];
sx q[0];
rz(-2.8594236) q[0];
sx q[0];
rz(-2.1532568) q[0];
rz(-pi) q[1];
rz(0.81014388) q[2];
sx q[2];
rz(-0.64254566) q[2];
sx q[2];
rz(2.0755529) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.12441382) q[1];
sx q[1];
rz(-0.22089566) q[1];
sx q[1];
rz(3.0934889) q[1];
rz(-pi) q[2];
rz(0.76830058) q[3];
sx q[3];
rz(-2.6242206) q[3];
sx q[3];
rz(-2.9381772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.540655) q[2];
sx q[2];
rz(-1.6161852) q[2];
sx q[2];
rz(1.8589004) q[2];
rz(0.013785275) q[3];
sx q[3];
rz(-2.0163592) q[3];
sx q[3];
rz(-1.3330601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29671058) q[0];
sx q[0];
rz(-1.7322576) q[0];
sx q[0];
rz(-2.5175293) q[0];
rz(-2.2804672) q[1];
sx q[1];
rz(-0.54665165) q[1];
sx q[1];
rz(0.13269592) q[1];
rz(-2.843905) q[2];
sx q[2];
rz(-0.33391446) q[2];
sx q[2];
rz(-3.1261974) q[2];
rz(0.91681963) q[3];
sx q[3];
rz(-2.0967757) q[3];
sx q[3];
rz(1.7158333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

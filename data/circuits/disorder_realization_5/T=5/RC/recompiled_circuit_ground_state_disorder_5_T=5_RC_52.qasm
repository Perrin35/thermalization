OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(-1.0125546) q[0];
sx q[0];
rz(-2.2192686) q[0];
rz(2.2946279) q[1];
sx q[1];
rz(-1.4743409) q[1];
sx q[1];
rz(-0.21811952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1058278) q[0];
sx q[0];
rz(-2.0313861) q[0];
sx q[0];
rz(-1.6732107) q[0];
x q[1];
rz(3.0384205) q[2];
sx q[2];
rz(-0.37984797) q[2];
sx q[2];
rz(-1.0641629) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74789819) q[1];
sx q[1];
rz(-1.5285349) q[1];
sx q[1];
rz(1.3512004) q[1];
rz(-pi) q[2];
rz(-2.0006337) q[3];
sx q[3];
rz(-1.6419193) q[3];
sx q[3];
rz(-1.5971368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.073079022) q[2];
sx q[2];
rz(-1.928669) q[2];
sx q[2];
rz(-0.53654137) q[2];
rz(1.0088629) q[3];
sx q[3];
rz(-0.77102414) q[3];
sx q[3];
rz(-2.3276276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.63758481) q[0];
sx q[0];
rz(-1.8300087) q[0];
sx q[0];
rz(0.017039321) q[0];
rz(1.0631961) q[1];
sx q[1];
rz(-0.76779643) q[1];
sx q[1];
rz(2.1868736) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9794036) q[0];
sx q[0];
rz(-1.6260176) q[0];
sx q[0];
rz(2.9121141) q[0];
rz(2.9770697) q[2];
sx q[2];
rz(-1.8662819) q[2];
sx q[2];
rz(0.72876677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3054637) q[1];
sx q[1];
rz(-1.8082128) q[1];
sx q[1];
rz(-1.8029455) q[1];
rz(-pi) q[2];
rz(1.9185136) q[3];
sx q[3];
rz(-3.0849815) q[3];
sx q[3];
rz(-2.268102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9118328) q[2];
sx q[2];
rz(-2.2288897) q[2];
sx q[2];
rz(2.0274053) q[2];
rz(1.1344502) q[3];
sx q[3];
rz(-2.9454234) q[3];
sx q[3];
rz(-2.9451059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.590362) q[0];
sx q[0];
rz(-2.1901665) q[0];
sx q[0];
rz(0.1828585) q[0];
rz(3.0602449) q[1];
sx q[1];
rz(-0.75129879) q[1];
sx q[1];
rz(-2.523211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8326022) q[0];
sx q[0];
rz(-1.6438577) q[0];
sx q[0];
rz(-2.354524) q[0];
rz(-pi) q[1];
rz(-0.19962991) q[2];
sx q[2];
rz(-2.3702894) q[2];
sx q[2];
rz(0.16801258) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32633852) q[1];
sx q[1];
rz(-0.36838057) q[1];
sx q[1];
rz(-1.9659564) q[1];
rz(2.2719273) q[3];
sx q[3];
rz(-1.0683064) q[3];
sx q[3];
rz(1.089407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1387834) q[2];
sx q[2];
rz(-1.323779) q[2];
sx q[2];
rz(2.7749824) q[2];
rz(1.2159411) q[3];
sx q[3];
rz(-1.1953019) q[3];
sx q[3];
rz(2.478157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4902896) q[0];
sx q[0];
rz(-1.2325352) q[0];
sx q[0];
rz(2.41462) q[0];
rz(-0.90826774) q[1];
sx q[1];
rz(-2.5636702) q[1];
sx q[1];
rz(2.0982826) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3573049) q[0];
sx q[0];
rz(-2.4803964) q[0];
sx q[0];
rz(-0.0088328514) q[0];
rz(-1.1953148) q[2];
sx q[2];
rz(-1.5133148) q[2];
sx q[2];
rz(2.6901746) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.336652) q[1];
sx q[1];
rz(-2.0139317) q[1];
sx q[1];
rz(-0.27733485) q[1];
rz(-2.1328385) q[3];
sx q[3];
rz(-1.2444454) q[3];
sx q[3];
rz(0.98451738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.72994453) q[2];
sx q[2];
rz(-2.8204155) q[2];
sx q[2];
rz(-1.3058861) q[2];
rz(-3.0236687) q[3];
sx q[3];
rz(-0.4685466) q[3];
sx q[3];
rz(-1.0606631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0333198) q[0];
sx q[0];
rz(-0.98618996) q[0];
sx q[0];
rz(-1.6075851) q[0];
rz(0.29574212) q[1];
sx q[1];
rz(-1.60138) q[1];
sx q[1];
rz(-1.6827513) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6687377) q[0];
sx q[0];
rz(-2.3928527) q[0];
sx q[0];
rz(1.4677901) q[0];
rz(2.2590911) q[2];
sx q[2];
rz(-0.46445981) q[2];
sx q[2];
rz(1.8732173) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0572059) q[1];
sx q[1];
rz(-0.73250341) q[1];
sx q[1];
rz(-1.2744348) q[1];
rz(-pi) q[2];
rz(-0.98990654) q[3];
sx q[3];
rz(-1.7647997) q[3];
sx q[3];
rz(-1.7542183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38833388) q[2];
sx q[2];
rz(-2.7497141) q[2];
sx q[2];
rz(0.64588109) q[2];
rz(-1.9258457) q[3];
sx q[3];
rz(-1.4589717) q[3];
sx q[3];
rz(1.3418044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011768613) q[0];
sx q[0];
rz(-0.83368603) q[0];
sx q[0];
rz(0.60761333) q[0];
rz(1.9629078) q[1];
sx q[1];
rz(-1.8557529) q[1];
sx q[1];
rz(0.36373055) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0640489) q[0];
sx q[0];
rz(-1.0564959) q[0];
sx q[0];
rz(0.27405996) q[0];
rz(0.046437101) q[2];
sx q[2];
rz(-1.2383467) q[2];
sx q[2];
rz(2.9118371) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1235001) q[1];
sx q[1];
rz(-0.82595968) q[1];
sx q[1];
rz(-0.1634181) q[1];
rz(-0.97447864) q[3];
sx q[3];
rz(-1.5913561) q[3];
sx q[3];
rz(2.1316776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58934775) q[2];
sx q[2];
rz(-0.10422464) q[2];
sx q[2];
rz(1.1927401) q[2];
rz(0.73181152) q[3];
sx q[3];
rz(-1.7164427) q[3];
sx q[3];
rz(1.4881136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28855395) q[0];
sx q[0];
rz(-0.1703425) q[0];
sx q[0];
rz(-1.5227675) q[0];
rz(1.6743926) q[1];
sx q[1];
rz(-1.3102945) q[1];
sx q[1];
rz(-0.67749643) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8554133) q[0];
sx q[0];
rz(-2.2660846) q[0];
sx q[0];
rz(-3.0762818) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67393731) q[2];
sx q[2];
rz(-2.8455685) q[2];
sx q[2];
rz(2.7715671) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.575857) q[1];
sx q[1];
rz(-1.9565554) q[1];
sx q[1];
rz(2.0292474) q[1];
x q[2];
rz(1.7301637) q[3];
sx q[3];
rz(-2.2969118) q[3];
sx q[3];
rz(2.1350525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7447394) q[2];
sx q[2];
rz(-2.2262959) q[2];
sx q[2];
rz(1.3057365) q[2];
rz(2.8030677) q[3];
sx q[3];
rz(-2.0207696) q[3];
sx q[3];
rz(-0.6849851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5386388) q[0];
sx q[0];
rz(-2.9260577) q[0];
sx q[0];
rz(-0.18390528) q[0];
rz(-1.6869847) q[1];
sx q[1];
rz(-0.55211663) q[1];
sx q[1];
rz(-2.6457381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496807) q[0];
sx q[0];
rz(-0.56654585) q[0];
sx q[0];
rz(-2.0828155) q[0];
rz(0.83424904) q[2];
sx q[2];
rz(-0.96820346) q[2];
sx q[2];
rz(-1.1863697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.715818) q[1];
sx q[1];
rz(-1.1905021) q[1];
sx q[1];
rz(0.38399314) q[1];
rz(-pi) q[2];
rz(3.0227376) q[3];
sx q[3];
rz(-1.0215852) q[3];
sx q[3];
rz(-1.0224179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.054606525) q[2];
sx q[2];
rz(-0.84422529) q[2];
sx q[2];
rz(-1.1620713) q[2];
rz(1.6880796) q[3];
sx q[3];
rz(-0.96265692) q[3];
sx q[3];
rz(1.8614205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42089713) q[0];
sx q[0];
rz(-1.1524042) q[0];
sx q[0];
rz(-2.5757117) q[0];
rz(0.3903009) q[1];
sx q[1];
rz(-1.5719599) q[1];
sx q[1];
rz(1.8338667) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6193251) q[0];
sx q[0];
rz(-0.59460708) q[0];
sx q[0];
rz(0.99647605) q[0];
x q[1];
rz(2.0686223) q[2];
sx q[2];
rz(-0.42440571) q[2];
sx q[2];
rz(-2.5967732) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3257287) q[1];
sx q[1];
rz(-2.1727409) q[1];
sx q[1];
rz(-1.8900327) q[1];
rz(1.7203749) q[3];
sx q[3];
rz(-1.6892551) q[3];
sx q[3];
rz(1.7427904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2719443) q[2];
sx q[2];
rz(-0.41676909) q[2];
sx q[2];
rz(1.0375674) q[2];
rz(-1.8218254) q[3];
sx q[3];
rz(-2.3128553) q[3];
sx q[3];
rz(2.493609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8808402) q[0];
sx q[0];
rz(-0.97761959) q[0];
sx q[0];
rz(-0.64224893) q[0];
rz(1.3308659) q[1];
sx q[1];
rz(-0.86527491) q[1];
sx q[1];
rz(2.9790402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35643043) q[0];
sx q[0];
rz(-0.9097865) q[0];
sx q[0];
rz(-2.6399773) q[0];
x q[1];
rz(0.75679512) q[2];
sx q[2];
rz(-2.6399603) q[2];
sx q[2];
rz(-2.121814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.049920883) q[1];
sx q[1];
rz(-2.0494048) q[1];
sx q[1];
rz(-0.19772999) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58329669) q[3];
sx q[3];
rz(-1.5044893) q[3];
sx q[3];
rz(-2.5401859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6315397) q[2];
sx q[2];
rz(-1.2756462) q[2];
sx q[2];
rz(-1.0408545) q[2];
rz(-2.1736274) q[3];
sx q[3];
rz(-2.0716045) q[3];
sx q[3];
rz(0.66106838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.80617245) q[0];
sx q[0];
rz(-1.5652884) q[0];
sx q[0];
rz(-0.096927222) q[0];
rz(-2.9087635) q[1];
sx q[1];
rz(-1.0284582) q[1];
sx q[1];
rz(-3.1340541) q[1];
rz(2.3940968) q[2];
sx q[2];
rz(-1.9608254) q[2];
sx q[2];
rz(-1.335782) q[2];
rz(2.6951058) q[3];
sx q[3];
rz(-0.75110186) q[3];
sx q[3];
rz(2.5121381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

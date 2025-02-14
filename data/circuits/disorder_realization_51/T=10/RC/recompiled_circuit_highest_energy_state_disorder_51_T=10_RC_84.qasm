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
rz(-2.812204) q[0];
sx q[0];
rz(-0.63151276) q[0];
sx q[0];
rz(3.1407177) q[0];
rz(2.5007091) q[1];
sx q[1];
rz(-2.1407318) q[1];
sx q[1];
rz(-0.34520087) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54291081) q[0];
sx q[0];
rz(-3.0474159) q[0];
sx q[0];
rz(-1.3046632) q[0];
rz(-pi) q[1];
rz(-2.7852374) q[2];
sx q[2];
rz(-2.7389067) q[2];
sx q[2];
rz(2.8860983) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4507323) q[1];
sx q[1];
rz(-1.1530515) q[1];
sx q[1];
rz(2.2295206) q[1];
x q[2];
rz(1.3561068) q[3];
sx q[3];
rz(-3.0654538) q[3];
sx q[3];
rz(-0.74387459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9125646) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(2.933617) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(-2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308554) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(0.99288565) q[0];
rz(1.317124) q[1];
sx q[1];
rz(-0.31610745) q[1];
sx q[1];
rz(2.3670926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3178783) q[0];
sx q[0];
rz(-1.5817989) q[0];
sx q[0];
rz(1.7490134) q[0];
rz(-pi) q[1];
rz(2.6235145) q[2];
sx q[2];
rz(-2.2642914) q[2];
sx q[2];
rz(-2.3665646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44821139) q[1];
sx q[1];
rz(-0.44084545) q[1];
sx q[1];
rz(0.54852672) q[1];
rz(-pi) q[2];
rz(-0.22529545) q[3];
sx q[3];
rz(-1.1114128) q[3];
sx q[3];
rz(0.53250203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2449067) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(2.1265105) q[2];
rz(2.8924938) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(-2.7740313) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86431137) q[0];
sx q[0];
rz(-2.4503777) q[0];
sx q[0];
rz(0.15750289) q[0];
rz(-1.0109673) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(0.13883042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15181825) q[0];
sx q[0];
rz(-0.4954557) q[0];
sx q[0];
rz(2.302343) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1558258) q[2];
sx q[2];
rz(-2.2425644) q[2];
sx q[2];
rz(2.6443554) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5032387) q[1];
sx q[1];
rz(-1.8625096) q[1];
sx q[1];
rz(-2.9364763) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5957727) q[3];
sx q[3];
rz(-2.8132317) q[3];
sx q[3];
rz(0.92121802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6802754) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(2.7834564) q[2];
rz(-0.67874587) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(1.0391191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045227483) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(-2.8367693) q[0];
rz(-0.40799704) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-0.92794424) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.613425) q[0];
sx q[0];
rz(-1.4425462) q[0];
sx q[0];
rz(-1.7665461) q[0];
x q[1];
rz(2.5408542) q[2];
sx q[2];
rz(-1.033412) q[2];
sx q[2];
rz(1.5058277) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.18670652) q[1];
sx q[1];
rz(-2.2518603) q[1];
sx q[1];
rz(1.7658556) q[1];
rz(2.0870861) q[3];
sx q[3];
rz(-2.1355503) q[3];
sx q[3];
rz(-1.3206583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9503595) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(2.9534269) q[2];
rz(-1.2763216) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(3.1222043) q[0];
rz(1.060932) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(1.2394989) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4864522) q[0];
sx q[0];
rz(-2.1391272) q[0];
sx q[0];
rz(0.33354946) q[0];
rz(2.4024348) q[2];
sx q[2];
rz(-1.0598759) q[2];
sx q[2];
rz(0.3581274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8722788) q[1];
sx q[1];
rz(-1.8511103) q[1];
sx q[1];
rz(-0.53086908) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89770384) q[3];
sx q[3];
rz(-0.60808676) q[3];
sx q[3];
rz(-1.2230108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1218607) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(1.3266374) q[2];
rz(-2.895368) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6061848) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(0.73915172) q[0];
rz(0.34573653) q[1];
sx q[1];
rz(-1.812499) q[1];
sx q[1];
rz(-2.2703222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10340583) q[0];
sx q[0];
rz(-1.6400918) q[0];
sx q[0];
rz(2.359576) q[0];
rz(-pi) q[1];
rz(-0.99314697) q[2];
sx q[2];
rz(-0.9075853) q[2];
sx q[2];
rz(-3.1301067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4193486) q[1];
sx q[1];
rz(-0.39545317) q[1];
sx q[1];
rz(1.3802746) q[1];
rz(-2.0454117) q[3];
sx q[3];
rz(-1.0132917) q[3];
sx q[3];
rz(1.6826009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31200108) q[2];
sx q[2];
rz(-0.63903725) q[2];
sx q[2];
rz(-0.87257067) q[2];
rz(1.568659) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(-2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0503814) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(0.76164371) q[0];
rz(-0.42539445) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(-2.2147307) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0935912) q[0];
sx q[0];
rz(-0.33080745) q[0];
sx q[0];
rz(-0.76865102) q[0];
rz(2.2144187) q[2];
sx q[2];
rz(-2.1774051) q[2];
sx q[2];
rz(1.7546897) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4372448) q[1];
sx q[1];
rz(-1.3489375) q[1];
sx q[1];
rz(-1.941844) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6377444) q[3];
sx q[3];
rz(-0.62426011) q[3];
sx q[3];
rz(-0.74893307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1130134) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(-2.2704303) q[2];
rz(0.29843676) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(-1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262064) q[0];
sx q[0];
rz(-0.50006777) q[0];
sx q[0];
rz(0.4183847) q[0];
rz(0.2977953) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(-0.13430886) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1537841) q[0];
sx q[0];
rz(-1.0414818) q[0];
sx q[0];
rz(-2.3321926) q[0];
x q[1];
rz(1.9761132) q[2];
sx q[2];
rz(-1.7421054) q[2];
sx q[2];
rz(-1.5754981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1301535) q[1];
sx q[1];
rz(-1.4362037) q[1];
sx q[1];
rz(2.9573739) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0487324) q[3];
sx q[3];
rz(-0.52937859) q[3];
sx q[3];
rz(-0.10806187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(-0.64201391) q[2];
rz(2.0783453) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(2.1003907) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6006271) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(0.11216057) q[0];
rz(1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(-2.1122011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4698668) q[0];
sx q[0];
rz(-0.70505667) q[0];
sx q[0];
rz(1.5211578) q[0];
x q[1];
rz(-0.12597398) q[2];
sx q[2];
rz(-0.306189) q[2];
sx q[2];
rz(2.9842215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16949108) q[1];
sx q[1];
rz(-1.6147366) q[1];
sx q[1];
rz(-1.8010587) q[1];
rz(-1.4791895) q[3];
sx q[3];
rz(-1.1044958) q[3];
sx q[3];
rz(-2.8331851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7182497) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(0.038012803) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(-3.0003701) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24899471) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(0.41879642) q[0];
rz(-1.8839802) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(0.14258252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2217956) q[0];
sx q[0];
rz(-2.9807973) q[0];
sx q[0];
rz(0.17531403) q[0];
rz(-pi) q[1];
rz(2.5129287) q[2];
sx q[2];
rz(-1.7970048) q[2];
sx q[2];
rz(-2.5369801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5113509) q[1];
sx q[1];
rz(-0.6410743) q[1];
sx q[1];
rz(-1.1180693) q[1];
rz(-pi) q[2];
rz(2.7975818) q[3];
sx q[3];
rz(-2.0075402) q[3];
sx q[3];
rz(-0.94319447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3451781) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(-0.26819116) q[2];
rz(1.7372519) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(-1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0161229) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(0.12915962) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(-1.2186981) q[2];
sx q[2];
rz(-2.0205971) q[2];
sx q[2];
rz(-1.481075) q[2];
rz(-0.066349647) q[3];
sx q[3];
rz(-1.8649615) q[3];
sx q[3];
rz(-2.1129114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

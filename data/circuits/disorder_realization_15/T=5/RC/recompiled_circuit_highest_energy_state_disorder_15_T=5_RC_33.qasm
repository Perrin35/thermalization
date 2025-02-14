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
rz(1.6871356) q[0];
sx q[0];
rz(-1.1633101) q[0];
sx q[0];
rz(-1.2679201) q[0];
rz(1.6549702) q[1];
sx q[1];
rz(-1.854719) q[1];
sx q[1];
rz(0.16965228) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171863) q[0];
sx q[0];
rz(-0.61530441) q[0];
sx q[0];
rz(-2.2236373) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0526161) q[2];
sx q[2];
rz(-1.3498326) q[2];
sx q[2];
rz(-2.882304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3006966) q[1];
sx q[1];
rz(-1.54956) q[1];
sx q[1];
rz(-1.9663215) q[1];
x q[2];
rz(1.2581014) q[3];
sx q[3];
rz(-1.7732244) q[3];
sx q[3];
rz(2.4656338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9156645) q[2];
sx q[2];
rz(-1.6253086) q[2];
sx q[2];
rz(3.1212659) q[2];
rz(2.9018371) q[3];
sx q[3];
rz(-2.8862947) q[3];
sx q[3];
rz(-2.3026626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13424419) q[0];
sx q[0];
rz(-0.61779314) q[0];
sx q[0];
rz(-2.0452621) q[0];
rz(1.4758551) q[1];
sx q[1];
rz(-1.9225537) q[1];
sx q[1];
rz(0.79442564) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61928328) q[0];
sx q[0];
rz(-1.2566152) q[0];
sx q[0];
rz(2.7854325) q[0];
rz(-pi) q[1];
rz(2.3874823) q[2];
sx q[2];
rz(-1.1271141) q[2];
sx q[2];
rz(-0.59690969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8960171) q[1];
sx q[1];
rz(-1.4652677) q[1];
sx q[1];
rz(0.42550931) q[1];
x q[2];
rz(2.3268324) q[3];
sx q[3];
rz(-1.4346806) q[3];
sx q[3];
rz(-0.55908113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.085792556) q[2];
sx q[2];
rz(-1.5726568) q[2];
sx q[2];
rz(0.22573486) q[2];
rz(-2.8977532) q[3];
sx q[3];
rz(-2.1832681) q[3];
sx q[3];
rz(-2.9624511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2589515) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(0.42695811) q[0];
rz(-2.7985165) q[1];
sx q[1];
rz(-1.1767358) q[1];
sx q[1];
rz(-2.8899736) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93676567) q[0];
sx q[0];
rz(-2.511095) q[0];
sx q[0];
rz(0.38602324) q[0];
rz(-pi) q[1];
rz(-2.919056) q[2];
sx q[2];
rz(-2.3250569) q[2];
sx q[2];
rz(-0.41513069) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3494663) q[1];
sx q[1];
rz(-1.9053962) q[1];
sx q[1];
rz(-1.3931331) q[1];
rz(-pi) q[2];
rz(1.6289457) q[3];
sx q[3];
rz(-2.4302501) q[3];
sx q[3];
rz(0.99866801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.130927) q[2];
sx q[2];
rz(-1.5074573) q[2];
sx q[2];
rz(0.44535401) q[2];
rz(-1.2238067) q[3];
sx q[3];
rz(-1.8755951) q[3];
sx q[3];
rz(-0.38770097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7447516) q[0];
sx q[0];
rz(-0.80533177) q[0];
sx q[0];
rz(-1.9824363) q[0];
rz(-1.2027488) q[1];
sx q[1];
rz(-1.9614599) q[1];
sx q[1];
rz(2.4045827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096367717) q[0];
sx q[0];
rz(-2.8857433) q[0];
sx q[0];
rz(-0.25080983) q[0];
x q[1];
rz(2.2463655) q[2];
sx q[2];
rz(-1.8631808) q[2];
sx q[2];
rz(1.043373) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9056978) q[1];
sx q[1];
rz(-1.7445843) q[1];
sx q[1];
rz(-0.64069616) q[1];
rz(2.8094711) q[3];
sx q[3];
rz(-0.17617036) q[3];
sx q[3];
rz(0.15819269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1272588) q[2];
sx q[2];
rz(-2.2517683) q[2];
sx q[2];
rz(0.8117525) q[2];
rz(0.74770606) q[3];
sx q[3];
rz(-0.85993189) q[3];
sx q[3];
rz(-3.0809793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6478445) q[0];
sx q[0];
rz(-1.0890549) q[0];
sx q[0];
rz(-1.3519979) q[0];
rz(2.3494675) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(2.1314714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8787694) q[0];
sx q[0];
rz(-2.4973625) q[0];
sx q[0];
rz(-0.99272002) q[0];
x q[1];
rz(-0.3926362) q[2];
sx q[2];
rz(-1.5199888) q[2];
sx q[2];
rz(-1.9399373) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3492584) q[1];
sx q[1];
rz(-2.5790303) q[1];
sx q[1];
rz(-2.9922561) q[1];
rz(0.95577464) q[3];
sx q[3];
rz(-0.20032756) q[3];
sx q[3];
rz(-2.1708787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.036756214) q[2];
sx q[2];
rz(-2.3827621) q[2];
sx q[2];
rz(1.7581615) q[2];
rz(-3.1116327) q[3];
sx q[3];
rz(-2.1432917) q[3];
sx q[3];
rz(1.2768607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6303915) q[0];
sx q[0];
rz(-2.7250405) q[0];
sx q[0];
rz(3.0101486) q[0];
rz(-2.1427872) q[1];
sx q[1];
rz(-1.9824332) q[1];
sx q[1];
rz(-1.3712032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5906487) q[0];
sx q[0];
rz(-1.2281795) q[0];
sx q[0];
rz(2.7638373) q[0];
x q[1];
rz(0.68626257) q[2];
sx q[2];
rz(-1.7289844) q[2];
sx q[2];
rz(-2.7686518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9483395) q[1];
sx q[1];
rz(-2.8242556) q[1];
sx q[1];
rz(2.0977519) q[1];
x q[2];
rz(-3.1375299) q[3];
sx q[3];
rz(-0.71089478) q[3];
sx q[3];
rz(1.9181741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3580631) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(1.0895458) q[2];
rz(0.36695668) q[3];
sx q[3];
rz(-0.4042545) q[3];
sx q[3];
rz(0.75016108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51339665) q[0];
sx q[0];
rz(-0.80444002) q[0];
sx q[0];
rz(-0.87271571) q[0];
rz(1.5645507) q[1];
sx q[1];
rz(-1.6460452) q[1];
sx q[1];
rz(2.4094792) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66595378) q[0];
sx q[0];
rz(-1.3697962) q[0];
sx q[0];
rz(-2.4253286) q[0];
x q[1];
rz(1.4462412) q[2];
sx q[2];
rz(-0.8890748) q[2];
sx q[2];
rz(-2.9602697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70591761) q[1];
sx q[1];
rz(-2.6443548) q[1];
sx q[1];
rz(-1.9051208) q[1];
rz(-pi) q[2];
x q[2];
rz(0.096166178) q[3];
sx q[3];
rz(-1.3447666) q[3];
sx q[3];
rz(-1.9440252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6096036) q[2];
sx q[2];
rz(-1.0010109) q[2];
sx q[2];
rz(-0.58094376) q[2];
rz(-2.0161435) q[3];
sx q[3];
rz(-1.9476798) q[3];
sx q[3];
rz(-1.9862991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693212) q[0];
sx q[0];
rz(-0.10675616) q[0];
sx q[0];
rz(-0.17053764) q[0];
rz(-1.8960309) q[1];
sx q[1];
rz(-1.6163369) q[1];
sx q[1];
rz(0.68444288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4675198) q[0];
sx q[0];
rz(-2.0430123) q[0];
sx q[0];
rz(-1.945482) q[0];
rz(-pi) q[1];
rz(-0.48284097) q[2];
sx q[2];
rz(-1.6833268) q[2];
sx q[2];
rz(1.7140599) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0674921) q[1];
sx q[1];
rz(-1.3252186) q[1];
sx q[1];
rz(2.8771993) q[1];
rz(-pi) q[2];
rz(-0.6489469) q[3];
sx q[3];
rz(-1.5228378) q[3];
sx q[3];
rz(3.0779882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.30274063) q[2];
sx q[2];
rz(-1.3417696) q[2];
sx q[2];
rz(2.0212685) q[2];
rz(0.66425792) q[3];
sx q[3];
rz(-0.40226007) q[3];
sx q[3];
rz(0.50814381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6054194) q[0];
sx q[0];
rz(-2.7689731) q[0];
sx q[0];
rz(-2.5033409) q[0];
rz(3.1328746) q[1];
sx q[1];
rz(-1.1161085) q[1];
sx q[1];
rz(-0.4062103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0779408) q[0];
sx q[0];
rz(-1.6839875) q[0];
sx q[0];
rz(1.8454396) q[0];
rz(-1.4582602) q[2];
sx q[2];
rz(-0.43607682) q[2];
sx q[2];
rz(-1.9364408) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6510424) q[1];
sx q[1];
rz(-0.79572751) q[1];
sx q[1];
rz(2.9814238) q[1];
rz(-0.33274098) q[3];
sx q[3];
rz(-0.077692835) q[3];
sx q[3];
rz(1.2487703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.093988769) q[2];
sx q[2];
rz(-2.0145907) q[2];
sx q[2];
rz(-0.36063933) q[2];
rz(2.4141198) q[3];
sx q[3];
rz(-0.68220264) q[3];
sx q[3];
rz(-0.31807652) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3280846) q[0];
sx q[0];
rz(-2.1961975) q[0];
sx q[0];
rz(-0.16950053) q[0];
rz(-2.4318579) q[1];
sx q[1];
rz(-1.7252445) q[1];
sx q[1];
rz(-1.08606) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5616578) q[0];
sx q[0];
rz(-2.5626963) q[0];
sx q[0];
rz(-2.7282342) q[0];
x q[1];
rz(0.31137054) q[2];
sx q[2];
rz(-1.1527921) q[2];
sx q[2];
rz(0.091146221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6352313) q[1];
sx q[1];
rz(-1.6690134) q[1];
sx q[1];
rz(1.0088831) q[1];
rz(0.38137718) q[3];
sx q[3];
rz(-1.3998195) q[3];
sx q[3];
rz(-1.7538296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3544932) q[2];
sx q[2];
rz(-0.084986173) q[2];
sx q[2];
rz(-0.3698012) q[2];
rz(1.9434816) q[3];
sx q[3];
rz(-1.5503649) q[3];
sx q[3];
rz(0.91355356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030466) q[0];
sx q[0];
rz(-1.1644762) q[0];
sx q[0];
rz(-1.2559011) q[0];
rz(-1.0176324) q[1];
sx q[1];
rz(-1.7435278) q[1];
sx q[1];
rz(1.7494038) q[1];
rz(3.0690003) q[2];
sx q[2];
rz(-1.6410927) q[2];
sx q[2];
rz(-2.1319364) q[2];
rz(-2.9774278) q[3];
sx q[3];
rz(-1.4420684) q[3];
sx q[3];
rz(1.3185929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

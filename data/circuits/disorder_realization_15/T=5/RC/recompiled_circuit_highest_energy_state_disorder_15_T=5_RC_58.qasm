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
rz(1.8736725) q[0];
rz(-1.4866225) q[1];
sx q[1];
rz(-1.2868737) q[1];
sx q[1];
rz(2.9719404) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21714144) q[0];
sx q[0];
rz(-1.0944774) q[0];
sx q[0];
rz(-0.40556559) q[0];
x q[1];
rz(1.1450139) q[2];
sx q[2];
rz(-0.55934956) q[2];
sx q[2];
rz(0.9445136) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.7387661) q[1];
sx q[1];
rz(-1.1753653) q[1];
sx q[1];
rz(-3.1185802) q[1];
rz(-pi) q[2];
rz(-2.9291487) q[3];
sx q[3];
rz(-1.8768969) q[3];
sx q[3];
rz(0.82992947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9156645) q[2];
sx q[2];
rz(-1.516284) q[2];
sx q[2];
rz(-3.1212659) q[2];
rz(-0.23975553) q[3];
sx q[3];
rz(-2.8862947) q[3];
sx q[3];
rz(0.83893004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13424419) q[0];
sx q[0];
rz(-0.61779314) q[0];
sx q[0];
rz(2.0452621) q[0];
rz(-1.4758551) q[1];
sx q[1];
rz(-1.9225537) q[1];
sx q[1];
rz(2.347167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4969032) q[0];
sx q[0];
rz(-0.4704501) q[0];
sx q[0];
rz(-0.75017922) q[0];
rz(0.99296221) q[2];
sx q[2];
rz(-2.2374399) q[2];
sx q[2];
rz(1.3571598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8640722) q[1];
sx q[1];
rz(-1.1478067) q[1];
sx q[1];
rz(-1.4550243) q[1];
rz(-2.3268324) q[3];
sx q[3];
rz(-1.4346806) q[3];
sx q[3];
rz(0.55908113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0558001) q[2];
sx q[2];
rz(-1.5726568) q[2];
sx q[2];
rz(-2.9158578) q[2];
rz(2.8977532) q[3];
sx q[3];
rz(-0.95832458) q[3];
sx q[3];
rz(-2.9624511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2589515) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(-2.7146345) q[0];
rz(-0.34307617) q[1];
sx q[1];
rz(-1.9648569) q[1];
sx q[1];
rz(0.25161904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31684012) q[0];
sx q[0];
rz(-1.7946294) q[0];
sx q[0];
rz(0.59455183) q[0];
x q[1];
rz(1.8015091) q[2];
sx q[2];
rz(-0.78015155) q[2];
sx q[2];
rz(3.0456269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9906944) q[1];
sx q[1];
rz(-0.37726918) q[1];
sx q[1];
rz(-0.47028415) q[1];
rz(-pi) q[2];
rz(3.0915458) q[3];
sx q[3];
rz(-0.86090961) q[3];
sx q[3];
rz(2.0662226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0106657) q[2];
sx q[2];
rz(-1.6341354) q[2];
sx q[2];
rz(0.44535401) q[2];
rz(-1.917786) q[3];
sx q[3];
rz(-1.2659975) q[3];
sx q[3];
rz(-0.38770097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39684108) q[0];
sx q[0];
rz(-0.80533177) q[0];
sx q[0];
rz(1.9824363) q[0];
rz(-1.2027488) q[1];
sx q[1];
rz(-1.1801327) q[1];
sx q[1];
rz(-2.4045827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4241959) q[0];
sx q[0];
rz(-1.6336461) q[0];
sx q[0];
rz(2.8934188) q[0];
rz(2.2463655) q[2];
sx q[2];
rz(-1.8631808) q[2];
sx q[2];
rz(-2.0982197) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.678469) q[1];
sx q[1];
rz(-0.94128525) q[1];
sx q[1];
rz(1.3552107) q[1];
rz(0.16672552) q[3];
sx q[3];
rz(-1.5136216) q[3];
sx q[3];
rz(-2.056332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1272588) q[2];
sx q[2];
rz(-0.88982439) q[2];
sx q[2];
rz(-2.3298402) q[2];
rz(-0.74770606) q[3];
sx q[3];
rz(-0.85993189) q[3];
sx q[3];
rz(3.0809793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6478445) q[0];
sx q[0];
rz(-2.0525377) q[0];
sx q[0];
rz(1.3519979) q[0];
rz(-2.3494675) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(1.0101213) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5631752) q[0];
sx q[0];
rz(-1.0437328) q[0];
sx q[0];
rz(0.38946797) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0094658) q[2];
sx q[2];
rz(-2.7458522) q[2];
sx q[2];
rz(-0.49116116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6163075) q[1];
sx q[1];
rz(-2.1263564) q[1];
sx q[1];
rz(1.6643334) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0249686) q[3];
sx q[3];
rz(-1.4075446) q[3];
sx q[3];
rz(-1.5463055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1048364) q[2];
sx q[2];
rz(-2.3827621) q[2];
sx q[2];
rz(1.3834312) q[2];
rz(-3.1116327) q[3];
sx q[3];
rz(-0.99830097) q[3];
sx q[3];
rz(-1.2768607) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6303915) q[0];
sx q[0];
rz(-2.7250405) q[0];
sx q[0];
rz(-0.1314441) q[0];
rz(-0.99880544) q[1];
sx q[1];
rz(-1.1591594) q[1];
sx q[1];
rz(-1.3712032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11267647) q[0];
sx q[0];
rz(-1.925591) q[0];
sx q[0];
rz(-1.9372069) q[0];
rz(-pi) q[1];
rz(2.8949685) q[2];
sx q[2];
rz(-0.70136753) q[2];
sx q[2];
rz(-1.3878617) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2590717) q[1];
sx q[1];
rz(-1.7283727) q[1];
sx q[1];
rz(1.2941918) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1375299) q[3];
sx q[3];
rz(-2.4306979) q[3];
sx q[3];
rz(-1.2234185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3580631) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(-2.0520468) q[2];
rz(-2.774636) q[3];
sx q[3];
rz(-2.7373382) q[3];
sx q[3];
rz(-0.75016108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628196) q[0];
sx q[0];
rz(-0.80444002) q[0];
sx q[0];
rz(0.87271571) q[0];
rz(1.577042) q[1];
sx q[1];
rz(-1.6460452) q[1];
sx q[1];
rz(-2.4094792) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0769147) q[0];
sx q[0];
rz(-2.2696662) q[0];
sx q[0];
rz(-1.3069606) q[0];
rz(-pi) q[1];
rz(2.9896834) q[2];
sx q[2];
rz(-0.6912125) q[2];
sx q[2];
rz(2.76413) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70591761) q[1];
sx q[1];
rz(-2.6443548) q[1];
sx q[1];
rz(-1.2364718) q[1];
rz(1.9663345) q[3];
sx q[3];
rz(-0.24531454) q[3];
sx q[3];
rz(-0.79110629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6096036) q[2];
sx q[2];
rz(-2.1405818) q[2];
sx q[2];
rz(-0.58094376) q[2];
rz(-1.1254492) q[3];
sx q[3];
rz(-1.1939129) q[3];
sx q[3];
rz(1.1552936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.67227143) q[0];
sx q[0];
rz(-3.0348365) q[0];
sx q[0];
rz(2.971055) q[0];
rz(1.2455617) q[1];
sx q[1];
rz(-1.5252557) q[1];
sx q[1];
rz(-0.68444288) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3869707) q[0];
sx q[0];
rz(-0.59387654) q[0];
sx q[0];
rz(0.62174364) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9028394) q[2];
sx q[2];
rz(-0.49477067) q[2];
sx q[2];
rz(-0.35428167) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2350133) q[1];
sx q[1];
rz(-2.7827127) q[1];
sx q[1];
rz(-0.76452629) q[1];
x q[2];
rz(-0.6489469) q[3];
sx q[3];
rz(-1.5228378) q[3];
sx q[3];
rz(3.0779882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.838852) q[2];
sx q[2];
rz(-1.799823) q[2];
sx q[2];
rz(-2.0212685) q[2];
rz(0.66425792) q[3];
sx q[3];
rz(-0.40226007) q[3];
sx q[3];
rz(0.50814381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53617322) q[0];
sx q[0];
rz(-2.7689731) q[0];
sx q[0];
rz(-0.63825178) q[0];
rz(0.0087180184) q[1];
sx q[1];
rz(-1.1161085) q[1];
sx q[1];
rz(-2.7353824) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4753301) q[0];
sx q[0];
rz(-1.8436369) q[0];
sx q[0];
rz(-0.11755875) q[0];
rz(-pi) q[1];
rz(-1.4582602) q[2];
sx q[2];
rz(-0.43607682) q[2];
sx q[2];
rz(1.2051518) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49055028) q[1];
sx q[1];
rz(-2.3458651) q[1];
sx q[1];
rz(2.9814238) q[1];
x q[2];
rz(2.8088517) q[3];
sx q[3];
rz(-3.0638998) q[3];
sx q[3];
rz(1.8928224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0476039) q[2];
sx q[2];
rz(-1.1270019) q[2];
sx q[2];
rz(-2.7809533) q[2];
rz(0.7274729) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(-0.31807652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3280846) q[0];
sx q[0];
rz(-0.94539517) q[0];
sx q[0];
rz(0.16950053) q[0];
rz(-2.4318579) q[1];
sx q[1];
rz(-1.4163481) q[1];
sx q[1];
rz(-2.0555326) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.062604) q[0];
sx q[0];
rz(-2.0955756) q[0];
sx q[0];
rz(1.3140509) q[0];
rz(-0.31137054) q[2];
sx q[2];
rz(-1.1527921) q[2];
sx q[2];
rz(-0.091146221) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12609657) q[1];
sx q[1];
rz(-1.0119146) q[1];
sx q[1];
rz(0.11591594) q[1];
x q[2];
rz(2.7602155) q[3];
sx q[3];
rz(-1.7417731) q[3];
sx q[3];
rz(-1.7538296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3544932) q[2];
sx q[2];
rz(-0.084986173) q[2];
sx q[2];
rz(-0.3698012) q[2];
rz(-1.1981111) q[3];
sx q[3];
rz(-1.5503649) q[3];
sx q[3];
rz(0.91355356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13854606) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(1.0176324) q[1];
sx q[1];
rz(-1.3980649) q[1];
sx q[1];
rz(-1.3921888) q[1];
rz(0.77059435) q[2];
sx q[2];
rz(-0.10100766) q[2];
sx q[2];
rz(-2.934692) q[2];
rz(-0.66987517) q[3];
sx q[3];
rz(-2.9333339) q[3];
sx q[3];
rz(0.40706709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

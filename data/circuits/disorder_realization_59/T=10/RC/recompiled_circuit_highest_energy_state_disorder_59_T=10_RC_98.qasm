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
rz(1.264313) q[0];
sx q[0];
rz(-2.8180583) q[0];
sx q[0];
rz(0.4488332) q[0];
rz(1.9857061) q[1];
sx q[1];
rz(-1.356025) q[1];
sx q[1];
rz(1.1745656) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6960775) q[0];
sx q[0];
rz(-0.68176523) q[0];
sx q[0];
rz(0.56607874) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5023793) q[2];
sx q[2];
rz(-1.2478154) q[2];
sx q[2];
rz(-3.0634865) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.60845637) q[1];
sx q[1];
rz(-1.605833) q[1];
sx q[1];
rz(2.1916981) q[1];
rz(0.51382084) q[3];
sx q[3];
rz(-1.4830477) q[3];
sx q[3];
rz(2.3140571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.590098) q[2];
sx q[2];
rz(-1.8415201) q[2];
sx q[2];
rz(-0.91280118) q[2];
rz(1.8707976) q[3];
sx q[3];
rz(-2.1015034) q[3];
sx q[3];
rz(-2.2916268) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0047282334) q[0];
sx q[0];
rz(-2.5623463) q[0];
sx q[0];
rz(2.7194523) q[0];
rz(-0.36960754) q[1];
sx q[1];
rz(-1.0969176) q[1];
sx q[1];
rz(0.94013989) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5135315) q[0];
sx q[0];
rz(-1.2413176) q[0];
sx q[0];
rz(0.0026436289) q[0];
rz(-0.11905382) q[2];
sx q[2];
rz(-1.040689) q[2];
sx q[2];
rz(-1.334545) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.861189) q[1];
sx q[1];
rz(-0.94377667) q[1];
sx q[1];
rz(0.95071937) q[1];
rz(-2.2475082) q[3];
sx q[3];
rz(-1.7853123) q[3];
sx q[3];
rz(-1.0971951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9422354) q[2];
sx q[2];
rz(-0.67023674) q[2];
sx q[2];
rz(2.2352236) q[2];
rz(-0.97936112) q[3];
sx q[3];
rz(-2.7370079) q[3];
sx q[3];
rz(-1.8766859) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449428) q[0];
sx q[0];
rz(-3.0746089) q[0];
sx q[0];
rz(-1.2445194) q[0];
rz(-2.7527346) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(-0.32274524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8334482) q[0];
sx q[0];
rz(-0.9128154) q[0];
sx q[0];
rz(1.4288155) q[0];
rz(-1.0110485) q[2];
sx q[2];
rz(-1.8343535) q[2];
sx q[2];
rz(-1.4562869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14798627) q[1];
sx q[1];
rz(-1.3393511) q[1];
sx q[1];
rz(2.3416421) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0067919) q[3];
sx q[3];
rz(-0.35405891) q[3];
sx q[3];
rz(-1.1969138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9129703) q[2];
sx q[2];
rz(-1.3801489) q[2];
sx q[2];
rz(1.7886394) q[2];
rz(-2.8690423) q[3];
sx q[3];
rz(-1.291178) q[3];
sx q[3];
rz(-1.4636309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54101855) q[0];
sx q[0];
rz(-0.32222846) q[0];
sx q[0];
rz(-1.898265) q[0];
rz(-2.2728032) q[1];
sx q[1];
rz(-0.96005762) q[1];
sx q[1];
rz(1.0702081) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7420827) q[0];
sx q[0];
rz(-2.1993981) q[0];
sx q[0];
rz(-0.6357155) q[0];
rz(2.3467873) q[2];
sx q[2];
rz(-1.2441976) q[2];
sx q[2];
rz(-2.2692533) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29954545) q[1];
sx q[1];
rz(-1.4562325) q[1];
sx q[1];
rz(1.24694) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34262212) q[3];
sx q[3];
rz(-1.3903769) q[3];
sx q[3];
rz(2.8460549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.59226817) q[2];
sx q[2];
rz(-0.78846875) q[2];
sx q[2];
rz(0.92174021) q[2];
rz(-2.8978469) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(2.8488979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9859966) q[0];
sx q[0];
rz(-1.8159741) q[0];
sx q[0];
rz(2.6407114) q[0];
rz(-2.0241375) q[1];
sx q[1];
rz(-1.3331579) q[1];
sx q[1];
rz(2.8270328) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524744) q[0];
sx q[0];
rz(-0.50273147) q[0];
sx q[0];
rz(-2.377654) q[0];
rz(-pi) q[1];
rz(0.0044745046) q[2];
sx q[2];
rz(-1.6142077) q[2];
sx q[2];
rz(-2.2501066) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.52175888) q[1];
sx q[1];
rz(-1.9653826) q[1];
sx q[1];
rz(2.16741) q[1];
rz(-pi) q[2];
rz(-1.8338449) q[3];
sx q[3];
rz(-0.42908731) q[3];
sx q[3];
rz(1.3671631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.306281) q[2];
sx q[2];
rz(-1.1113144) q[2];
sx q[2];
rz(1.8446946) q[2];
rz(-2.6563307) q[3];
sx q[3];
rz(-1.5596215) q[3];
sx q[3];
rz(-1.1135134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49419633) q[0];
sx q[0];
rz(-1.2911456) q[0];
sx q[0];
rz(3.0418292) q[0];
rz(-1.2707204) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(1.7521923) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6806599) q[0];
sx q[0];
rz(-1.7353199) q[0];
sx q[0];
rz(2.8712832) q[0];
x q[1];
rz(0.54047853) q[2];
sx q[2];
rz(-0.20217824) q[2];
sx q[2];
rz(2.0433189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2459316) q[1];
sx q[1];
rz(-0.97890039) q[1];
sx q[1];
rz(0.81879692) q[1];
rz(-pi) q[2];
rz(-0.78310599) q[3];
sx q[3];
rz(-1.4811188) q[3];
sx q[3];
rz(0.41148666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6483267) q[2];
sx q[2];
rz(-2.6706225) q[2];
sx q[2];
rz(-2.2512482) q[2];
rz(-1.9503615) q[3];
sx q[3];
rz(-1.8343364) q[3];
sx q[3];
rz(2.2395535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7118199) q[0];
sx q[0];
rz(-0.044782488) q[0];
sx q[0];
rz(3.1088767) q[0];
rz(0.50321594) q[1];
sx q[1];
rz(-1.0992173) q[1];
sx q[1];
rz(-3.018766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8211201) q[0];
sx q[0];
rz(-2.3373465) q[0];
sx q[0];
rz(-0.45566092) q[0];
x q[1];
rz(-0.57374222) q[2];
sx q[2];
rz(-2.5649568) q[2];
sx q[2];
rz(-1.5337616) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6325153) q[1];
sx q[1];
rz(-1.4715403) q[1];
sx q[1];
rz(-2.0429913) q[1];
rz(-0.35251725) q[3];
sx q[3];
rz(-1.4488701) q[3];
sx q[3];
rz(2.8059949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34933019) q[2];
sx q[2];
rz(-2.1749039) q[2];
sx q[2];
rz(1.1518504) q[2];
rz(2.2564015) q[3];
sx q[3];
rz(-1.7722292) q[3];
sx q[3];
rz(1.7159897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952878) q[0];
sx q[0];
rz(-2.2080244) q[0];
sx q[0];
rz(0.92920148) q[0];
rz(-0.62905351) q[1];
sx q[1];
rz(-1.355143) q[1];
sx q[1];
rz(-2.3984875) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5261054) q[0];
sx q[0];
rz(-1.1243404) q[0];
sx q[0];
rz(-0.087421405) q[0];
rz(-pi) q[1];
rz(-0.488398) q[2];
sx q[2];
rz(-1.6553859) q[2];
sx q[2];
rz(1.1957439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23452246) q[1];
sx q[1];
rz(-1.5046232) q[1];
sx q[1];
rz(-0.41445865) q[1];
rz(-pi) q[2];
rz(2.1984897) q[3];
sx q[3];
rz(-2.4022958) q[3];
sx q[3];
rz(3.0969248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5762081) q[2];
sx q[2];
rz(-2.3614063) q[2];
sx q[2];
rz(-0.72594491) q[2];
rz(2.7706326) q[3];
sx q[3];
rz(-1.5434664) q[3];
sx q[3];
rz(0.36271873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2711656) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(-0.56394947) q[0];
rz(-1.9922527) q[1];
sx q[1];
rz(-2.1795858) q[1];
sx q[1];
rz(-0.74251485) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71096984) q[0];
sx q[0];
rz(-1.1046243) q[0];
sx q[0];
rz(-0.52180334) q[0];
rz(-0.41461035) q[2];
sx q[2];
rz(-1.3002204) q[2];
sx q[2];
rz(2.0913578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2143605) q[1];
sx q[1];
rz(-0.4402658) q[1];
sx q[1];
rz(-3.0378067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56470676) q[3];
sx q[3];
rz(-2.2032524) q[3];
sx q[3];
rz(2.6868827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5138862) q[2];
sx q[2];
rz(-1.0739001) q[2];
sx q[2];
rz(-2.7970496) q[2];
rz(-2.6978317) q[3];
sx q[3];
rz(-2.0659476) q[3];
sx q[3];
rz(-1.3226604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76170707) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(0.66201061) q[0];
rz(-0.33540353) q[1];
sx q[1];
rz(-1.4619724) q[1];
sx q[1];
rz(0.9368771) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10231384) q[0];
sx q[0];
rz(-2.1697308) q[0];
sx q[0];
rz(2.7014669) q[0];
x q[1];
rz(-0.023597864) q[2];
sx q[2];
rz(-1.840971) q[2];
sx q[2];
rz(-0.23575704) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.42439991) q[1];
sx q[1];
rz(-0.48312995) q[1];
sx q[1];
rz(2.4880954) q[1];
x q[2];
rz(-2.2912737) q[3];
sx q[3];
rz(-0.26278824) q[3];
sx q[3];
rz(2.9966054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4235437) q[2];
sx q[2];
rz(-0.23423883) q[2];
sx q[2];
rz(-0.49918175) q[2];
rz(1.4623803) q[3];
sx q[3];
rz(-1.1464109) q[3];
sx q[3];
rz(1.4596938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.725631) q[0];
sx q[0];
rz(-1.9209296) q[0];
sx q[0];
rz(-0.85335535) q[0];
rz(0.60494963) q[1];
sx q[1];
rz(-1.9974983) q[1];
sx q[1];
rz(-2.6430184) q[1];
rz(-2.5685979) q[2];
sx q[2];
rz(-1.4481432) q[2];
sx q[2];
rz(-1.0021423) q[2];
rz(1.3868757) q[3];
sx q[3];
rz(-2.0652323) q[3];
sx q[3];
rz(-0.064814322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(0.19658495) q[0];
sx q[0];
rz(-2.6994446) q[0];
sx q[0];
rz(-0.86139876) q[0];
rz(-1.6317033) q[1];
sx q[1];
rz(-1.3246526) q[1];
sx q[1];
rz(-3.1060001) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1382768) q[0];
sx q[0];
rz(-0.39130032) q[0];
sx q[0];
rz(2.6075493) q[0];
rz(-1.6794559) q[2];
sx q[2];
rz(-1.8755967) q[2];
sx q[2];
rz(2.5632312) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.030499996) q[1];
sx q[1];
rz(-0.72734088) q[1];
sx q[1];
rz(3.0877581) q[1];
rz(-pi) q[2];
rz(1.54856) q[3];
sx q[3];
rz(-2.4614868) q[3];
sx q[3];
rz(0.82397991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55277905) q[2];
sx q[2];
rz(-1.3014883) q[2];
sx q[2];
rz(-1.0944875) q[2];
rz(1.6257446) q[3];
sx q[3];
rz(-0.87259126) q[3];
sx q[3];
rz(-2.998013) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0482408) q[0];
sx q[0];
rz(-2.1513262) q[0];
sx q[0];
rz(0.37407237) q[0];
rz(0.85809842) q[1];
sx q[1];
rz(-0.94081196) q[1];
sx q[1];
rz(-2.0657952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8825622) q[0];
sx q[0];
rz(-2.9963065) q[0];
sx q[0];
rz(2.3502716) q[0];
rz(-pi) q[1];
rz(1.3316989) q[2];
sx q[2];
rz(-2.1054724) q[2];
sx q[2];
rz(-2.3102674) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3871284) q[1];
sx q[1];
rz(-0.417388) q[1];
sx q[1];
rz(2.8512958) q[1];
x q[2];
rz(0.76794736) q[3];
sx q[3];
rz(-2.30184) q[3];
sx q[3];
rz(2.6382382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6263803) q[2];
sx q[2];
rz(-2.7084646) q[2];
sx q[2];
rz(-1.0832146) q[2];
rz(1.2927239) q[3];
sx q[3];
rz(-1.1336528) q[3];
sx q[3];
rz(-2.2129272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2556297) q[0];
sx q[0];
rz(-0.77115458) q[0];
sx q[0];
rz(-2.6166925) q[0];
rz(-1.4595754) q[1];
sx q[1];
rz(-2.4039905) q[1];
sx q[1];
rz(2.9958013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12851579) q[0];
sx q[0];
rz(-1.4798963) q[0];
sx q[0];
rz(-0.10581067) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1728233) q[2];
sx q[2];
rz(-1.1175691) q[2];
sx q[2];
rz(-2.9211958) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1555712) q[1];
sx q[1];
rz(-1.1386445) q[1];
sx q[1];
rz(-1.6701103) q[1];
rz(-2.1946241) q[3];
sx q[3];
rz(-0.83443975) q[3];
sx q[3];
rz(-2.8744227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8968481) q[2];
sx q[2];
rz(-1.9795828) q[2];
sx q[2];
rz(-0.16540089) q[2];
rz(2.2798955) q[3];
sx q[3];
rz(-0.69178897) q[3];
sx q[3];
rz(-2.944788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2794613) q[0];
sx q[0];
rz(-2.617351) q[0];
sx q[0];
rz(-2.4133546) q[0];
rz(-0.32314745) q[1];
sx q[1];
rz(-1.0341945) q[1];
sx q[1];
rz(-1.039215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6704039) q[0];
sx q[0];
rz(-1.9946596) q[0];
sx q[0];
rz(0.43760646) q[0];
rz(-pi) q[1];
rz(-0.28113725) q[2];
sx q[2];
rz(-2.1875811) q[2];
sx q[2];
rz(-2.5532818) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0919623) q[1];
sx q[1];
rz(-2.6481682) q[1];
sx q[1];
rz(2.9454524) q[1];
rz(-pi) q[2];
rz(2.3023241) q[3];
sx q[3];
rz(-2.2429562) q[3];
sx q[3];
rz(1.2109001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49179658) q[2];
sx q[2];
rz(-1.0504664) q[2];
sx q[2];
rz(0.51898471) q[2];
rz(1.0745878) q[3];
sx q[3];
rz(-1.855775) q[3];
sx q[3];
rz(-0.49526596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.0786667) q[0];
sx q[0];
rz(-0.92200297) q[0];
sx q[0];
rz(-2.9700188) q[0];
rz(-1.0942787) q[1];
sx q[1];
rz(-1.3010052) q[1];
sx q[1];
rz(-0.23460728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6162524) q[0];
sx q[0];
rz(-1.8484637) q[0];
sx q[0];
rz(-0.13292952) q[0];
x q[1];
rz(-0.50723664) q[2];
sx q[2];
rz(-1.5248858) q[2];
sx q[2];
rz(-3.1243589) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68069862) q[1];
sx q[1];
rz(-2.1354224) q[1];
sx q[1];
rz(-2.2959397) q[1];
rz(2.4314546) q[3];
sx q[3];
rz(-1.2810859) q[3];
sx q[3];
rz(0.75677815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91308633) q[2];
sx q[2];
rz(-1.9501016) q[2];
sx q[2];
rz(1.7388434) q[2];
rz(-0.62471041) q[3];
sx q[3];
rz(-2.1731845) q[3];
sx q[3];
rz(-1.7984084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45873555) q[0];
sx q[0];
rz(-0.0077489297) q[0];
sx q[0];
rz(0.32753456) q[0];
rz(0.028118357) q[1];
sx q[1];
rz(-1.1437806) q[1];
sx q[1];
rz(-0.99172529) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2153127) q[0];
sx q[0];
rz(-0.87020391) q[0];
sx q[0];
rz(0.69045569) q[0];
x q[1];
rz(2.9155511) q[2];
sx q[2];
rz(-1.6767354) q[2];
sx q[2];
rz(0.49913479) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8918482) q[1];
sx q[1];
rz(-1.6127024) q[1];
sx q[1];
rz(0.81467198) q[1];
rz(-pi) q[2];
x q[2];
rz(3.087849) q[3];
sx q[3];
rz(-0.47131594) q[3];
sx q[3];
rz(1.1137247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2345387) q[2];
sx q[2];
rz(-1.3776642) q[2];
sx q[2];
rz(-1.482359) q[2];
rz(-1.0659069) q[3];
sx q[3];
rz(-1.1803455) q[3];
sx q[3];
rz(1.0970241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3858353) q[0];
sx q[0];
rz(-0.11180728) q[0];
sx q[0];
rz(2.8461611) q[0];
rz(-0.23751986) q[1];
sx q[1];
rz(-0.93370456) q[1];
sx q[1];
rz(-2.3942153) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45493653) q[0];
sx q[0];
rz(-1.0191518) q[0];
sx q[0];
rz(-2.174189) q[0];
x q[1];
rz(2.105999) q[2];
sx q[2];
rz(-1.0135302) q[2];
sx q[2];
rz(-2.0095428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.805757) q[1];
sx q[1];
rz(-1.079396) q[1];
sx q[1];
rz(2.7193428) q[1];
x q[2];
rz(-1.9456057) q[3];
sx q[3];
rz(-2.0829933) q[3];
sx q[3];
rz(0.95694816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6682917) q[2];
sx q[2];
rz(-1.95582) q[2];
sx q[2];
rz(0.2549003) q[2];
rz(0.21236803) q[3];
sx q[3];
rz(-2.2771211) q[3];
sx q[3];
rz(-3.1406241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2905228) q[0];
sx q[0];
rz(-1.2238598) q[0];
sx q[0];
rz(0.12390027) q[0];
rz(2.2663785) q[1];
sx q[1];
rz(-2.3085935) q[1];
sx q[1];
rz(-0.55878729) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4717455) q[0];
sx q[0];
rz(-0.59006834) q[0];
sx q[0];
rz(1.1099932) q[0];
rz(0.28382878) q[2];
sx q[2];
rz(-2.0528194) q[2];
sx q[2];
rz(1.1555156) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9001101) q[1];
sx q[1];
rz(-2.0962976) q[1];
sx q[1];
rz(1.8837758) q[1];
x q[2];
rz(2.2053917) q[3];
sx q[3];
rz(-1.9125328) q[3];
sx q[3];
rz(0.5638916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0957886) q[2];
sx q[2];
rz(-0.49171058) q[2];
sx q[2];
rz(2.4208505) q[2];
rz(-2.7785684) q[3];
sx q[3];
rz(-1.9178773) q[3];
sx q[3];
rz(-1.5117517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95661288) q[0];
sx q[0];
rz(-2.9473801) q[0];
sx q[0];
rz(0.089056253) q[0];
rz(0.36674276) q[1];
sx q[1];
rz(-2.2373503) q[1];
sx q[1];
rz(-1.4595703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37292591) q[0];
sx q[0];
rz(-0.67250508) q[0];
sx q[0];
rz(1.9422533) q[0];
x q[1];
rz(-2.4581681) q[2];
sx q[2];
rz(-2.0439842) q[2];
sx q[2];
rz(-3.0101208) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68946099) q[1];
sx q[1];
rz(-2.2575827) q[1];
sx q[1];
rz(2.6584714) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59883786) q[3];
sx q[3];
rz(-1.9396848) q[3];
sx q[3];
rz(-1.2303331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9350932) q[2];
sx q[2];
rz(-1.3672028) q[2];
sx q[2];
rz(1.8915141) q[2];
rz(-1.2262723) q[3];
sx q[3];
rz(-0.63392249) q[3];
sx q[3];
rz(0.63898501) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6298237) q[0];
sx q[0];
rz(-2.0084232) q[0];
sx q[0];
rz(-2.0015707) q[0];
rz(-2.5584768) q[1];
sx q[1];
rz(-1.6551599) q[1];
sx q[1];
rz(-0.66782943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7191737) q[0];
sx q[0];
rz(-1.2454709) q[0];
sx q[0];
rz(0.76119411) q[0];
x q[1];
rz(-2.9029151) q[2];
sx q[2];
rz(-2.2584923) q[2];
sx q[2];
rz(-2.4872045) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8357667) q[1];
sx q[1];
rz(-2.0025616) q[1];
sx q[1];
rz(-2.9271896) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6891081) q[3];
sx q[3];
rz(-0.84558949) q[3];
sx q[3];
rz(3.0304327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84393152) q[2];
sx q[2];
rz(-0.49489489) q[2];
sx q[2];
rz(0.75221357) q[2];
rz(-0.080605896) q[3];
sx q[3];
rz(-1.7919431) q[3];
sx q[3];
rz(0.54752553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053631) q[0];
sx q[0];
rz(-1.6148051) q[0];
sx q[0];
rz(-0.057407277) q[0];
rz(-2.2930131) q[1];
sx q[1];
rz(-0.72863693) q[1];
sx q[1];
rz(-1.4562664) q[1];
rz(0.66203881) q[2];
sx q[2];
rz(-1.0218191) q[2];
sx q[2];
rz(-1.54984) q[2];
rz(-1.2101444) q[3];
sx q[3];
rz(-1.1359884) q[3];
sx q[3];
rz(2.8154919) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

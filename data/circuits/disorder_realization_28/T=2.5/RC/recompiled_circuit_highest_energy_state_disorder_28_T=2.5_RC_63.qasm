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
rz(-2.3856491) q[0];
sx q[0];
rz(-0.33060253) q[0];
sx q[0];
rz(-2.8688353) q[0];
rz(2.7859712) q[1];
sx q[1];
rz(1.0300809) q[1];
sx q[1];
rz(7.8539943) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8518675) q[0];
sx q[0];
rz(-0.12061943) q[0];
sx q[0];
rz(-0.9319181) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5982389) q[2];
sx q[2];
rz(-2.5961868) q[2];
sx q[2];
rz(-2.0871833) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94704365) q[1];
sx q[1];
rz(-1.6131372) q[1];
sx q[1];
rz(-2.4621099) q[1];
rz(-pi) q[2];
rz(-2.8924283) q[3];
sx q[3];
rz(-2.7621885) q[3];
sx q[3];
rz(-0.5642341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.97751272) q[2];
sx q[2];
rz(-1.432632) q[2];
sx q[2];
rz(-2.254503) q[2];
rz(-1.8767493) q[3];
sx q[3];
rz(-1.0586459) q[3];
sx q[3];
rz(-3.1177055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3905268) q[0];
sx q[0];
rz(-2.2845415) q[0];
sx q[0];
rz(2.8476025) q[0];
rz(1.7901621) q[1];
sx q[1];
rz(-1.0963115) q[1];
sx q[1];
rz(1.9757804) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0675177) q[0];
sx q[0];
rz(-0.39108983) q[0];
sx q[0];
rz(2.142201) q[0];
x q[1];
rz(-1.2109257) q[2];
sx q[2];
rz(-3.0635298) q[2];
sx q[2];
rz(-3.0541354) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68968523) q[1];
sx q[1];
rz(-2.4396585) q[1];
sx q[1];
rz(2.4268389) q[1];
rz(-0.43710093) q[3];
sx q[3];
rz(-0.61449285) q[3];
sx q[3];
rz(1.75625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68842781) q[2];
sx q[2];
rz(-0.63899779) q[2];
sx q[2];
rz(-1.2028018) q[2];
rz(1.5996251) q[3];
sx q[3];
rz(-1.8034233) q[3];
sx q[3];
rz(-0.28551027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0315392) q[0];
sx q[0];
rz(-0.078131214) q[0];
sx q[0];
rz(-1.1745289) q[0];
rz(-0.9749167) q[1];
sx q[1];
rz(-0.87923032) q[1];
sx q[1];
rz(2.6851795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.31118) q[0];
sx q[0];
rz(-1.5390656) q[0];
sx q[0];
rz(-1.8422442) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8936797) q[2];
sx q[2];
rz(-1.9447902) q[2];
sx q[2];
rz(-0.1318814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2172981) q[1];
sx q[1];
rz(-0.84554377) q[1];
sx q[1];
rz(2.3582703) q[1];
rz(-pi) q[2];
rz(0.47243709) q[3];
sx q[3];
rz(-2.2294913) q[3];
sx q[3];
rz(-2.6078826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6156893) q[2];
sx q[2];
rz(-0.38280767) q[2];
sx q[2];
rz(-0.75695401) q[2];
rz(-0.01903875) q[3];
sx q[3];
rz(-1.2913387) q[3];
sx q[3];
rz(-2.3058057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86209908) q[0];
sx q[0];
rz(-2.4308496) q[0];
sx q[0];
rz(-0.94974649) q[0];
rz(-2.9277335) q[1];
sx q[1];
rz(-0.93997926) q[1];
sx q[1];
rz(-2.5423539) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70940269) q[0];
sx q[0];
rz(-0.66973842) q[0];
sx q[0];
rz(-1.7959005) q[0];
rz(-1.2798115) q[2];
sx q[2];
rz(-2.6070234) q[2];
sx q[2];
rz(-0.9415516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1018925) q[1];
sx q[1];
rz(-1.3624316) q[1];
sx q[1];
rz(-1.031967) q[1];
rz(2.9327389) q[3];
sx q[3];
rz(-0.59126544) q[3];
sx q[3];
rz(-2.8610736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.63138258) q[2];
sx q[2];
rz(-1.6333132) q[2];
sx q[2];
rz(0.82376662) q[2];
rz(-1.0907762) q[3];
sx q[3];
rz(-0.80086509) q[3];
sx q[3];
rz(0.66528475) q[3];
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
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7432231) q[0];
sx q[0];
rz(-0.38341612) q[0];
sx q[0];
rz(-2.1552591) q[0];
rz(-2.8979454) q[1];
sx q[1];
rz(-1.2657093) q[1];
sx q[1];
rz(-0.15959127) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361317) q[0];
sx q[0];
rz(-1.2572932) q[0];
sx q[0];
rz(-2.3132597) q[0];
rz(-2.1417732) q[2];
sx q[2];
rz(-2.2717064) q[2];
sx q[2];
rz(-0.69055218) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21503285) q[1];
sx q[1];
rz(-1.6590282) q[1];
sx q[1];
rz(-2.0162321) q[1];
x q[2];
rz(2.2714204) q[3];
sx q[3];
rz(-1.7990094) q[3];
sx q[3];
rz(0.96503497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7470982) q[2];
sx q[2];
rz(-1.2622086) q[2];
sx q[2];
rz(-2.0060284) q[2];
rz(2.6731532) q[3];
sx q[3];
rz(-1.9197542) q[3];
sx q[3];
rz(-2.2170317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.121345) q[0];
sx q[0];
rz(-1.0806885) q[0];
sx q[0];
rz(-0.24170804) q[0];
rz(-1.7860335) q[1];
sx q[1];
rz(-0.87239289) q[1];
sx q[1];
rz(-2.2579069) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6970405) q[0];
sx q[0];
rz(-1.0844106) q[0];
sx q[0];
rz(2.4598511) q[0];
x q[1];
rz(1.8275798) q[2];
sx q[2];
rz(-1.2681539) q[2];
sx q[2];
rz(-3.1336477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7979398) q[1];
sx q[1];
rz(-2.6561497) q[1];
sx q[1];
rz(1.9992827) q[1];
rz(-pi) q[2];
rz(-3.1413636) q[3];
sx q[3];
rz(-1.159824) q[3];
sx q[3];
rz(2.4966893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9200865) q[2];
sx q[2];
rz(-0.93044996) q[2];
sx q[2];
rz(0.42631701) q[2];
rz(-0.91655556) q[3];
sx q[3];
rz(-1.5519658) q[3];
sx q[3];
rz(2.0385273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.351848) q[0];
sx q[0];
rz(-1.177657) q[0];
sx q[0];
rz(1.1258997) q[0];
rz(-3.0806091) q[1];
sx q[1];
rz(-2.1911502) q[1];
sx q[1];
rz(-2.1263863) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7506152) q[0];
sx q[0];
rz(-0.96834194) q[0];
sx q[0];
rz(1.0573122) q[0];
rz(1.3638968) q[2];
sx q[2];
rz(-2.7963429) q[2];
sx q[2];
rz(0.35790863) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0132844) q[1];
sx q[1];
rz(-2.7375712) q[1];
sx q[1];
rz(1.728251) q[1];
x q[2];
rz(3.0281248) q[3];
sx q[3];
rz(-1.5312636) q[3];
sx q[3];
rz(3.0276445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31275493) q[2];
sx q[2];
rz(-1.2042896) q[2];
sx q[2];
rz(0.22605669) q[2];
rz(2.1894646) q[3];
sx q[3];
rz(-3.003037) q[3];
sx q[3];
rz(0.8086732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1299745) q[0];
sx q[0];
rz(-1.8267153) q[0];
sx q[0];
rz(0.34039482) q[0];
rz(3.0420692) q[1];
sx q[1];
rz(-1.1025905) q[1];
sx q[1];
rz(-1.5955101) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9700546) q[0];
sx q[0];
rz(-2.4405789) q[0];
sx q[0];
rz(-0.51999493) q[0];
rz(-3.0853188) q[2];
sx q[2];
rz(-1.84308) q[2];
sx q[2];
rz(2.9416566) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9118285) q[1];
sx q[1];
rz(-0.10783261) q[1];
sx q[1];
rz(-2.7417438) q[1];
rz(2.9603954) q[3];
sx q[3];
rz(-0.46162039) q[3];
sx q[3];
rz(0.68440729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.978329) q[2];
sx q[2];
rz(-1.949387) q[2];
sx q[2];
rz(-2.1400129) q[2];
rz(-1.7892276) q[3];
sx q[3];
rz(-2.6632301) q[3];
sx q[3];
rz(-1.9549687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96420646) q[0];
sx q[0];
rz(-1.0564251) q[0];
sx q[0];
rz(-0.00038432234) q[0];
rz(2.7035825) q[1];
sx q[1];
rz(-1.6338467) q[1];
sx q[1];
rz(0.45723525) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5075574) q[0];
sx q[0];
rz(-2.9762161) q[0];
sx q[0];
rz(1.2654113) q[0];
rz(-2.1887378) q[2];
sx q[2];
rz(-2.0597201) q[2];
sx q[2];
rz(3.0143154) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36553069) q[1];
sx q[1];
rz(-1.9161822) q[1];
sx q[1];
rz(-1.8074492) q[1];
x q[2];
rz(2.6351686) q[3];
sx q[3];
rz(-1.2791787) q[3];
sx q[3];
rz(-0.17403655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75927258) q[2];
sx q[2];
rz(-0.98733941) q[2];
sx q[2];
rz(-1.6299204) q[2];
rz(-3.0509389) q[3];
sx q[3];
rz(-1.0414618) q[3];
sx q[3];
rz(2.4618728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35336211) q[0];
sx q[0];
rz(-1.3557949) q[0];
sx q[0];
rz(-0.28025383) q[0];
rz(1.7578112) q[1];
sx q[1];
rz(-2.6774112) q[1];
sx q[1];
rz(-1.9162477) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2458151) q[0];
sx q[0];
rz(-1.4207463) q[0];
sx q[0];
rz(-1.4729985) q[0];
rz(-pi) q[1];
rz(2.2544075) q[2];
sx q[2];
rz(-1.2204183) q[2];
sx q[2];
rz(-2.7840707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30825057) q[1];
sx q[1];
rz(-2.3099515) q[1];
sx q[1];
rz(-1.2396481) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2810542) q[3];
sx q[3];
rz(-1.2803161) q[3];
sx q[3];
rz(-1.2742467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9489991) q[2];
sx q[2];
rz(-1.5263564) q[2];
sx q[2];
rz(-2.4601649) q[2];
rz(2.0415908) q[3];
sx q[3];
rz(-2.2097578) q[3];
sx q[3];
rz(-1.8345691) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.451402) q[0];
sx q[0];
rz(-0.9342397) q[0];
sx q[0];
rz(-1.1697212) q[0];
rz(-2.7611217) q[1];
sx q[1];
rz(-2.4741551) q[1];
sx q[1];
rz(-2.9173775) q[1];
rz(-2.8773814) q[2];
sx q[2];
rz(-0.88709863) q[2];
sx q[2];
rz(-1.8567793) q[2];
rz(-0.43222618) q[3];
sx q[3];
rz(-1.0226915) q[3];
sx q[3];
rz(2.6998479) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

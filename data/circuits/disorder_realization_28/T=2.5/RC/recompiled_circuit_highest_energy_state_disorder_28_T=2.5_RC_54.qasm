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
rz(0.75594354) q[0];
sx q[0];
rz(3.4721952) q[0];
sx q[0];
rz(9.1520206) q[0];
rz(2.7859712) q[1];
sx q[1];
rz(1.0300809) q[1];
sx q[1];
rz(7.8539943) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2897252) q[0];
sx q[0];
rz(-0.12061943) q[0];
sx q[0];
rz(0.9319181) q[0];
x q[1];
rz(1.2414957) q[2];
sx q[2];
rz(-1.1277756) q[2];
sx q[2];
rz(-1.7274513) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.57140162) q[1];
sx q[1];
rz(-0.68059151) q[1];
sx q[1];
rz(0.067318214) q[1];
rz(-pi) q[2];
rz(1.4727888) q[3];
sx q[3];
rz(-1.2036754) q[3];
sx q[3];
rz(-2.844732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97751272) q[2];
sx q[2];
rz(-1.7089607) q[2];
sx q[2];
rz(-0.88708964) q[2];
rz(-1.2648434) q[3];
sx q[3];
rz(-2.0829468) q[3];
sx q[3];
rz(0.023887159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.3905268) q[0];
sx q[0];
rz(-0.85705119) q[0];
sx q[0];
rz(0.29399011) q[0];
rz(-1.7901621) q[1];
sx q[1];
rz(-2.0452812) q[1];
sx q[1];
rz(1.9757804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0675177) q[0];
sx q[0];
rz(-0.39108983) q[0];
sx q[0];
rz(-2.142201) q[0];
rz(-pi) q[1];
rz(-0.027539081) q[2];
sx q[2];
rz(-1.4977432) q[2];
sx q[2];
rz(2.6932584) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4519074) q[1];
sx q[1];
rz(-2.4396585) q[1];
sx q[1];
rz(2.4268389) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8610624) q[3];
sx q[3];
rz(-1.0212047) q[3];
sx q[3];
rz(0.86587469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.68842781) q[2];
sx q[2];
rz(-0.63899779) q[2];
sx q[2];
rz(1.2028018) q[2];
rz(-1.5996251) q[3];
sx q[3];
rz(-1.3381693) q[3];
sx q[3];
rz(2.8560824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11005345) q[0];
sx q[0];
rz(-3.0634614) q[0];
sx q[0];
rz(1.9670638) q[0];
rz(0.9749167) q[1];
sx q[1];
rz(-0.87923032) q[1];
sx q[1];
rz(0.45641315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62688359) q[0];
sx q[0];
rz(-2.868342) q[0];
sx q[0];
rz(1.4529626) q[0];
x q[1];
rz(1.012031) q[2];
sx q[2];
rz(-0.44545275) q[2];
sx q[2];
rz(-0.73778668) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3756936) q[1];
sx q[1];
rz(-1.0119034) q[1];
sx q[1];
rz(-2.2431348) q[1];
x q[2];
rz(-2.2863488) q[3];
sx q[3];
rz(-1.9389098) q[3];
sx q[3];
rz(2.4076622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5259033) q[2];
sx q[2];
rz(-0.38280767) q[2];
sx q[2];
rz(0.75695401) q[2];
rz(-0.01903875) q[3];
sx q[3];
rz(-1.8502539) q[3];
sx q[3];
rz(-0.835787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86209908) q[0];
sx q[0];
rz(-0.71074301) q[0];
sx q[0];
rz(-2.1918462) q[0];
rz(-0.21385916) q[1];
sx q[1];
rz(-2.2016134) q[1];
sx q[1];
rz(-2.5423539) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1480162) q[0];
sx q[0];
rz(-0.92088078) q[0];
sx q[0];
rz(2.9666569) q[0];
rz(-pi) q[1];
rz(-1.2798115) q[2];
sx q[2];
rz(-2.6070234) q[2];
sx q[2];
rz(-0.9415516) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1018925) q[1];
sx q[1];
rz(-1.779161) q[1];
sx q[1];
rz(-2.1096257) q[1];
x q[2];
rz(2.5604519) q[3];
sx q[3];
rz(-1.4549644) q[3];
sx q[3];
rz(2.0254962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5102101) q[2];
sx q[2];
rz(-1.5082794) q[2];
sx q[2];
rz(-2.317826) q[2];
rz(2.0508164) q[3];
sx q[3];
rz(-2.3407276) q[3];
sx q[3];
rz(-0.66528475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7432231) q[0];
sx q[0];
rz(-0.38341612) q[0];
sx q[0];
rz(-2.1552591) q[0];
rz(-0.24364722) q[1];
sx q[1];
rz(-1.8758834) q[1];
sx q[1];
rz(2.9820014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8411113) q[0];
sx q[0];
rz(-2.2694025) q[0];
sx q[0];
rz(-0.41450702) q[0];
x q[1];
rz(0.56964376) q[2];
sx q[2];
rz(-0.87235556) q[2];
sx q[2];
rz(-1.6676355) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9683384) q[1];
sx q[1];
rz(-2.6880777) q[1];
sx q[1];
rz(-1.7733002) q[1];
rz(-pi) q[2];
rz(-2.8466264) q[3];
sx q[3];
rz(-0.89184232) q[3];
sx q[3];
rz(2.7243638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7470982) q[2];
sx q[2];
rz(-1.8793841) q[2];
sx q[2];
rz(2.0060284) q[2];
rz(2.6731532) q[3];
sx q[3];
rz(-1.9197542) q[3];
sx q[3];
rz(-2.2170317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0202476) q[0];
sx q[0];
rz(-2.0609042) q[0];
sx q[0];
rz(-0.24170804) q[0];
rz(-1.3555591) q[1];
sx q[1];
rz(-2.2691998) q[1];
sx q[1];
rz(-2.2579069) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4445522) q[0];
sx q[0];
rz(-1.0844106) q[0];
sx q[0];
rz(-0.68174151) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3140129) q[2];
sx q[2];
rz(-1.8734387) q[2];
sx q[2];
rz(0.0079449991) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3436529) q[1];
sx q[1];
rz(-2.6561497) q[1];
sx q[1];
rz(-1.14231) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9817686) q[3];
sx q[3];
rz(-1.5710063) q[3];
sx q[3];
rz(-2.2157912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22150618) q[2];
sx q[2];
rz(-2.2111427) q[2];
sx q[2];
rz(-2.7152756) q[2];
rz(-0.91655556) q[3];
sx q[3];
rz(-1.5896268) q[3];
sx q[3];
rz(-2.0385273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78974462) q[0];
sx q[0];
rz(-1.9639356) q[0];
sx q[0];
rz(-2.015693) q[0];
rz(3.0806091) q[1];
sx q[1];
rz(-2.1911502) q[1];
sx q[1];
rz(-1.0152063) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7495817) q[0];
sx q[0];
rz(-0.77031743) q[0];
sx q[0];
rz(2.5213741) q[0];
x q[1];
rz(3.0678441) q[2];
sx q[2];
rz(-1.2332067) q[2];
sx q[2];
rz(-3.003157) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0132844) q[1];
sx q[1];
rz(-2.7375712) q[1];
sx q[1];
rz(1.728251) q[1];
x q[2];
rz(-3.0281248) q[3];
sx q[3];
rz(-1.610329) q[3];
sx q[3];
rz(-0.11394812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31275493) q[2];
sx q[2];
rz(-1.9373031) q[2];
sx q[2];
rz(-0.22605669) q[2];
rz(2.1894646) q[3];
sx q[3];
rz(-0.13855562) q[3];
sx q[3];
rz(2.3329195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1299745) q[0];
sx q[0];
rz(-1.8267153) q[0];
sx q[0];
rz(-0.34039482) q[0];
rz(-0.099523425) q[1];
sx q[1];
rz(-1.1025905) q[1];
sx q[1];
rz(1.5460825) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9700546) q[0];
sx q[0];
rz(-2.4405789) q[0];
sx q[0];
rz(0.51999493) q[0];
rz(-pi) q[1];
x q[1];
rz(1.298102) q[2];
sx q[2];
rz(-1.6249949) q[2];
sx q[2];
rz(1.3860089) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2297641) q[1];
sx q[1];
rz(-3.03376) q[1];
sx q[1];
rz(-2.7417438) q[1];
x q[2];
rz(0.45507064) q[3];
sx q[3];
rz(-1.4904456) q[3];
sx q[3];
rz(0.72380607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.16326363) q[2];
sx q[2];
rz(-1.1922057) q[2];
sx q[2];
rz(-1.0015798) q[2];
rz(-1.352365) q[3];
sx q[3];
rz(-2.6632301) q[3];
sx q[3];
rz(-1.1866239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773862) q[0];
sx q[0];
rz(-1.0564251) q[0];
sx q[0];
rz(0.00038432234) q[0];
rz(0.43801019) q[1];
sx q[1];
rz(-1.507746) q[1];
sx q[1];
rz(-2.6843574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9033636) q[0];
sx q[0];
rz(-1.6203124) q[0];
sx q[0];
rz(-1.4129421) q[0];
x q[1];
rz(0.82797956) q[2];
sx q[2];
rz(-0.767602) q[2];
sx q[2];
rz(1.1141675) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8548566) q[1];
sx q[1];
rz(-1.7932307) q[1];
sx q[1];
rz(0.35448928) q[1];
x q[2];
rz(1.9014539) q[3];
sx q[3];
rz(-2.0539502) q[3];
sx q[3];
rz(1.5866947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3823201) q[2];
sx q[2];
rz(-0.98733941) q[2];
sx q[2];
rz(-1.6299204) q[2];
rz(-0.090653732) q[3];
sx q[3];
rz(-1.0414618) q[3];
sx q[3];
rz(-2.4618728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35336211) q[0];
sx q[0];
rz(-1.7857977) q[0];
sx q[0];
rz(-0.28025383) q[0];
rz(1.7578112) q[1];
sx q[1];
rz(-2.6774112) q[1];
sx q[1];
rz(-1.9162477) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8266164) q[0];
sx q[0];
rz(-0.17890636) q[0];
sx q[0];
rz(-0.57344435) q[0];
rz(-0.44048788) q[2];
sx q[2];
rz(-0.9357399) q[2];
sx q[2];
rz(-1.4859413) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8333421) q[1];
sx q[1];
rz(-0.83164117) q[1];
sx q[1];
rz(-1.9019446) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30237301) q[3];
sx q[3];
rz(-1.8480715) q[3];
sx q[3];
rz(-0.38173166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9489991) q[2];
sx q[2];
rz(-1.5263564) q[2];
sx q[2];
rz(-0.68142778) q[2];
rz(1.1000018) q[3];
sx q[3];
rz(-0.93183485) q[3];
sx q[3];
rz(-1.8345691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.451402) q[0];
sx q[0];
rz(-2.207353) q[0];
sx q[0];
rz(1.9718715) q[0];
rz(-2.7611217) q[1];
sx q[1];
rz(-2.4741551) q[1];
sx q[1];
rz(-2.9173775) q[1];
rz(-0.86974551) q[2];
sx q[2];
rz(-1.7746584) q[2];
sx q[2];
rz(-0.45523582) q[2];
rz(-0.97888234) q[3];
sx q[3];
rz(-1.2051734) q[3];
sx q[3];
rz(-1.7766374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

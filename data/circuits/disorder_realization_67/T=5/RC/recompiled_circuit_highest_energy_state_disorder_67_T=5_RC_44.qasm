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
rz(-0.29210583) q[0];
sx q[0];
rz(-0.6120683) q[0];
sx q[0];
rz(0.050961343) q[0];
rz(1.2904957) q[1];
sx q[1];
rz(-1.2761071) q[1];
sx q[1];
rz(2.2955503) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1109306) q[0];
sx q[0];
rz(-1.7314859) q[0];
sx q[0];
rz(-0.35344214) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1943075) q[2];
sx q[2];
rz(-1.3007015) q[2];
sx q[2];
rz(-0.17284753) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1345299) q[1];
sx q[1];
rz(-1.0872528) q[1];
sx q[1];
rz(-1.8731432) q[1];
x q[2];
rz(-1.5791164) q[3];
sx q[3];
rz(-1.4589615) q[3];
sx q[3];
rz(0.40598265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35752615) q[2];
sx q[2];
rz(-2.1937723) q[2];
sx q[2];
rz(-0.86304647) q[2];
rz(-1.1132025) q[3];
sx q[3];
rz(-0.92156711) q[3];
sx q[3];
rz(-0.79871261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5384101) q[0];
sx q[0];
rz(-0.69284678) q[0];
sx q[0];
rz(2.8811654) q[0];
rz(0.32568359) q[1];
sx q[1];
rz(-0.98418701) q[1];
sx q[1];
rz(-2.8338285) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64806077) q[0];
sx q[0];
rz(-2.2770398) q[0];
sx q[0];
rz(0.34428455) q[0];
rz(-pi) q[1];
rz(-0.4701647) q[2];
sx q[2];
rz(-1.2855296) q[2];
sx q[2];
rz(-0.72848749) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5883397) q[1];
sx q[1];
rz(-0.25811895) q[1];
sx q[1];
rz(1.2023167) q[1];
rz(-pi) q[2];
x q[2];
rz(0.006853718) q[3];
sx q[3];
rz(-0.93539507) q[3];
sx q[3];
rz(1.8346661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8327568) q[2];
sx q[2];
rz(-2.2726161) q[2];
sx q[2];
rz(-0.55574179) q[2];
rz(-0.51221383) q[3];
sx q[3];
rz(-0.54849505) q[3];
sx q[3];
rz(0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.888805) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(0.74932253) q[0];
rz(-2.0454171) q[1];
sx q[1];
rz(-1.7671894) q[1];
sx q[1];
rz(-0.38415092) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2040981) q[0];
sx q[0];
rz(-0.9956313) q[0];
sx q[0];
rz(-1.8857486) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60127705) q[2];
sx q[2];
rz(-1.7814751) q[2];
sx q[2];
rz(2.7803068) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2652581) q[1];
sx q[1];
rz(-1.5175845) q[1];
sx q[1];
rz(1.268178) q[1];
rz(1.2189193) q[3];
sx q[3];
rz(-2.4259544) q[3];
sx q[3];
rz(-1.972769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7595547) q[2];
sx q[2];
rz(-1.3052772) q[2];
sx q[2];
rz(-1.6273512) q[2];
rz(-0.73005992) q[3];
sx q[3];
rz(-2.4000013) q[3];
sx q[3];
rz(-2.7470284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6877614) q[0];
sx q[0];
rz(-1.7173959) q[0];
sx q[0];
rz(3.0815531) q[0];
rz(2.2278348) q[1];
sx q[1];
rz(-0.91122952) q[1];
sx q[1];
rz(0.93889108) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1969827) q[0];
sx q[0];
rz(-1.6337578) q[0];
sx q[0];
rz(1.5374166) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.082415103) q[2];
sx q[2];
rz(-1.6019099) q[2];
sx q[2];
rz(1.4084202) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1891494) q[1];
sx q[1];
rz(-1.5464142) q[1];
sx q[1];
rz(-2.0588751) q[1];
rz(-pi) q[2];
rz(-2.3653829) q[3];
sx q[3];
rz(-0.84899711) q[3];
sx q[3];
rz(-1.2152023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.011220304) q[2];
sx q[2];
rz(-0.29272407) q[2];
sx q[2];
rz(-0.65227738) q[2];
rz(-1.8782714) q[3];
sx q[3];
rz(-1.5247366) q[3];
sx q[3];
rz(-1.9704845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071282722) q[0];
sx q[0];
rz(-1.902453) q[0];
sx q[0];
rz(2.5791445) q[0];
rz(1.3485472) q[1];
sx q[1];
rz(-2.8548073) q[1];
sx q[1];
rz(-0.58132344) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5347915) q[0];
sx q[0];
rz(-1.1282217) q[0];
sx q[0];
rz(-2.3032196) q[0];
rz(-1.5346709) q[2];
sx q[2];
rz(-1.110552) q[2];
sx q[2];
rz(-0.66260105) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5116252) q[1];
sx q[1];
rz(-1.9251578) q[1];
sx q[1];
rz(-2.0767861) q[1];
x q[2];
rz(1.6147037) q[3];
sx q[3];
rz(-2.2128519) q[3];
sx q[3];
rz(1.2894443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44744667) q[2];
sx q[2];
rz(-1.7116825) q[2];
sx q[2];
rz(2.4549129) q[2];
rz(-0.40003362) q[3];
sx q[3];
rz(-1.2587849) q[3];
sx q[3];
rz(0.011628477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7069063) q[0];
sx q[0];
rz(-2.1676368) q[0];
sx q[0];
rz(-2.6913225) q[0];
rz(-2.2606668) q[1];
sx q[1];
rz(-1.9073146) q[1];
sx q[1];
rz(-1.1164104) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0562949) q[0];
sx q[0];
rz(-1.3816815) q[0];
sx q[0];
rz(-2.3448639) q[0];
x q[1];
rz(-0.35043244) q[2];
sx q[2];
rz(-1.3964126) q[2];
sx q[2];
rz(-2.7728105) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9444869) q[1];
sx q[1];
rz(-2.53105) q[1];
sx q[1];
rz(-1.6050299) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0209337) q[3];
sx q[3];
rz(-1.880247) q[3];
sx q[3];
rz(-1.7904953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.474596) q[2];
sx q[2];
rz(-2.907739) q[2];
sx q[2];
rz(-2.4901566) q[2];
rz(2.3598119) q[3];
sx q[3];
rz(-1.4835446) q[3];
sx q[3];
rz(-0.93530161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5930138) q[0];
sx q[0];
rz(-0.23603708) q[0];
sx q[0];
rz(-2.6759942) q[0];
rz(-1.9904526) q[1];
sx q[1];
rz(-1.8955889) q[1];
sx q[1];
rz(-1.8021072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27522408) q[0];
sx q[0];
rz(-1.0491519) q[0];
sx q[0];
rz(2.0417968) q[0];
rz(-pi) q[1];
rz(0.8546245) q[2];
sx q[2];
rz(-3.0833088) q[2];
sx q[2];
rz(0.56266498) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7581343) q[1];
sx q[1];
rz(-1.7900568) q[1];
sx q[1];
rz(-3.0644528) q[1];
rz(-pi) q[2];
rz(2.2364113) q[3];
sx q[3];
rz(-1.2318805) q[3];
sx q[3];
rz(-1.9012332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0307978) q[2];
sx q[2];
rz(-1.3862627) q[2];
sx q[2];
rz(0.01290713) q[2];
rz(3.0068093) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(-0.24136647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6808692) q[0];
sx q[0];
rz(-3.0496106) q[0];
sx q[0];
rz(-1.9665834) q[0];
rz(0.75042945) q[1];
sx q[1];
rz(-1.3576327) q[1];
sx q[1];
rz(-2.7480385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9703044) q[0];
sx q[0];
rz(-1.7455202) q[0];
sx q[0];
rz(-2.1615684) q[0];
rz(-pi) q[1];
rz(-0.15109328) q[2];
sx q[2];
rz(-2.4799621) q[2];
sx q[2];
rz(-0.90666319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9266415) q[1];
sx q[1];
rz(-1.1424173) q[1];
sx q[1];
rz(3.0207816) q[1];
rz(-pi) q[2];
rz(2.7435947) q[3];
sx q[3];
rz(-1.0915712) q[3];
sx q[3];
rz(-2.3204539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63341081) q[2];
sx q[2];
rz(-2.0076624) q[2];
sx q[2];
rz(-2.9122162) q[2];
rz(-0.091382787) q[3];
sx q[3];
rz(-1.6978076) q[3];
sx q[3];
rz(0.55575931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2185739) q[0];
sx q[0];
rz(-0.78892437) q[0];
sx q[0];
rz(0.58513418) q[0];
rz(-1.8642037) q[1];
sx q[1];
rz(-1.9265415) q[1];
sx q[1];
rz(-2.9871984) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0994249) q[0];
sx q[0];
rz(-1.6800361) q[0];
sx q[0];
rz(-0.43719284) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4677708) q[2];
sx q[2];
rz(-0.5267064) q[2];
sx q[2];
rz(0.4835085) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2626434) q[1];
sx q[1];
rz(-0.41987881) q[1];
sx q[1];
rz(0.68563571) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2326689) q[3];
sx q[3];
rz(-0.49175374) q[3];
sx q[3];
rz(1.7724747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1207017) q[2];
sx q[2];
rz(-0.88940826) q[2];
sx q[2];
rz(-0.54023877) q[2];
rz(2.7464416) q[3];
sx q[3];
rz(-0.65427798) q[3];
sx q[3];
rz(-1.7925709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2111135) q[0];
sx q[0];
rz(-1.8607288) q[0];
sx q[0];
rz(1.6092009) q[0];
rz(-2.2156175) q[1];
sx q[1];
rz(-1.4264868) q[1];
sx q[1];
rz(-1.4769953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4187936) q[0];
sx q[0];
rz(-0.72615756) q[0];
sx q[0];
rz(3.1405057) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8468644) q[2];
sx q[2];
rz(-1.2090313) q[2];
sx q[2];
rz(-2.8579762) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.6546911) q[1];
sx q[1];
rz(-2.4729252) q[1];
sx q[1];
rz(-2.1441441) q[1];
x q[2];
rz(2.6981101) q[3];
sx q[3];
rz(-2.1135619) q[3];
sx q[3];
rz(-0.51668985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8605139) q[2];
sx q[2];
rz(-1.3382567) q[2];
sx q[2];
rz(0.64858428) q[2];
rz(-0.72566882) q[3];
sx q[3];
rz(-0.23894335) q[3];
sx q[3];
rz(1.749595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3758748) q[0];
sx q[0];
rz(-2.1633042) q[0];
sx q[0];
rz(-2.3194585) q[0];
rz(-1.0646461) q[1];
sx q[1];
rz(-1.4551661) q[1];
sx q[1];
rz(-1.0600769) q[1];
rz(-2.9072472) q[2];
sx q[2];
rz(-0.6498944) q[2];
sx q[2];
rz(1.9882974) q[2];
rz(2.3271046) q[3];
sx q[3];
rz(-1.1753488) q[3];
sx q[3];
rz(-1.4594441) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

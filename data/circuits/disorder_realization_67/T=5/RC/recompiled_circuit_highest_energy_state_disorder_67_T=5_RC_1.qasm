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
rz(2.8494868) q[0];
sx q[0];
rz(3.7536609) q[0];
sx q[0];
rz(6.232224) q[0];
rz(-1.851097) q[1];
sx q[1];
rz(-1.8654856) q[1];
sx q[1];
rz(-2.2955503) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.030662) q[0];
sx q[0];
rz(-1.4101068) q[0];
sx q[0];
rz(-2.7881505) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2162391) q[2];
sx q[2];
rz(-0.45956372) q[2];
sx q[2];
rz(1.1499576) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1345299) q[1];
sx q[1];
rz(-1.0872528) q[1];
sx q[1];
rz(1.2684494) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11183864) q[3];
sx q[3];
rz(-1.5625283) q[3];
sx q[3];
rz(-1.9777075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7840665) q[2];
sx q[2];
rz(-2.1937723) q[2];
sx q[2];
rz(2.2785462) q[2];
rz(-2.0283902) q[3];
sx q[3];
rz(-0.92156711) q[3];
sx q[3];
rz(-2.34288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6031826) q[0];
sx q[0];
rz(-0.69284678) q[0];
sx q[0];
rz(0.26042724) q[0];
rz(2.8159091) q[1];
sx q[1];
rz(-0.98418701) q[1];
sx q[1];
rz(-0.30776417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69410283) q[0];
sx q[0];
rz(-1.3110975) q[0];
sx q[0];
rz(-0.83456852) q[0];
rz(-0.57450104) q[2];
sx q[2];
rz(-0.5443474) q[2];
sx q[2];
rz(-2.8050204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1733117) q[1];
sx q[1];
rz(-1.3303583) q[1];
sx q[1];
rz(3.0467826) q[1];
x q[2];
rz(0.006853718) q[3];
sx q[3];
rz(-0.93539507) q[3];
sx q[3];
rz(-1.3069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8327568) q[2];
sx q[2];
rz(-2.2726161) q[2];
sx q[2];
rz(2.5858509) q[2];
rz(-0.51221383) q[3];
sx q[3];
rz(-2.5930976) q[3];
sx q[3];
rz(-0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.25278768) q[0];
sx q[0];
rz(-2.7894809) q[0];
sx q[0];
rz(0.74932253) q[0];
rz(2.0454171) q[1];
sx q[1];
rz(-1.3744033) q[1];
sx q[1];
rz(-0.38415092) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2040981) q[0];
sx q[0];
rz(-2.1459614) q[0];
sx q[0];
rz(1.8857486) q[0];
x q[1];
rz(-0.36142771) q[2];
sx q[2];
rz(-2.5088032) q[2];
sx q[2];
rz(-0.91362) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2652581) q[1];
sx q[1];
rz(-1.5175845) q[1];
sx q[1];
rz(-1.268178) q[1];
rz(-pi) q[2];
rz(-1.2189193) q[3];
sx q[3];
rz(-0.71563827) q[3];
sx q[3];
rz(1.1688237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7595547) q[2];
sx q[2];
rz(-1.3052772) q[2];
sx q[2];
rz(-1.5142415) q[2];
rz(-0.73005992) q[3];
sx q[3];
rz(-0.74159139) q[3];
sx q[3];
rz(-0.39456427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45383129) q[0];
sx q[0];
rz(-1.4241968) q[0];
sx q[0];
rz(0.060039595) q[0];
rz(-0.91375786) q[1];
sx q[1];
rz(-0.91122952) q[1];
sx q[1];
rz(0.93889108) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37591463) q[0];
sx q[0];
rz(-1.6041099) q[0];
sx q[0];
rz(-0.062996431) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0591776) q[2];
sx q[2];
rz(-1.6019099) q[2];
sx q[2];
rz(1.7331725) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36870391) q[1];
sx q[1];
rz(-2.0587173) q[1];
sx q[1];
rz(0.027603961) q[1];
rz(0.77620971) q[3];
sx q[3];
rz(-0.84899711) q[3];
sx q[3];
rz(1.9263903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.011220304) q[2];
sx q[2];
rz(-2.8488686) q[2];
sx q[2];
rz(-0.65227738) q[2];
rz(1.2633213) q[3];
sx q[3];
rz(-1.5247366) q[3];
sx q[3];
rz(-1.9704845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0703099) q[0];
sx q[0];
rz(-1.2391397) q[0];
sx q[0];
rz(-0.5624482) q[0];
rz(-1.3485472) q[1];
sx q[1];
rz(-2.8548073) q[1];
sx q[1];
rz(0.58132344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59636694) q[0];
sx q[0];
rz(-2.2196182) q[0];
sx q[0];
rz(-2.5741386) q[0];
rz(-pi) q[1];
rz(-1.5346709) q[2];
sx q[2];
rz(-2.0310406) q[2];
sx q[2];
rz(-2.4789916) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.86921924) q[1];
sx q[1];
rz(-1.0989184) q[1];
sx q[1];
rz(2.7414338) q[1];
x q[2];
rz(3.082959) q[3];
sx q[3];
rz(-2.4982493) q[3];
sx q[3];
rz(-1.7789121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.694146) q[2];
sx q[2];
rz(-1.4299102) q[2];
sx q[2];
rz(2.4549129) q[2];
rz(2.741559) q[3];
sx q[3];
rz(-1.2587849) q[3];
sx q[3];
rz(0.011628477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7069063) q[0];
sx q[0];
rz(-0.97395581) q[0];
sx q[0];
rz(-2.6913225) q[0];
rz(0.88092583) q[1];
sx q[1];
rz(-1.2342781) q[1];
sx q[1];
rz(1.1164104) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4661144) q[0];
sx q[0];
rz(-0.7921392) q[0];
sx q[0];
rz(-1.303543) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47412761) q[2];
sx q[2];
rz(-0.38981405) q[2];
sx q[2];
rz(-1.6451943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3456381) q[1];
sx q[1];
rz(-1.5511723) q[1];
sx q[1];
rz(0.96052891) q[1];
rz(-pi) q[2];
rz(2.0209337) q[3];
sx q[3];
rz(-1.2613457) q[3];
sx q[3];
rz(1.3510973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.474596) q[2];
sx q[2];
rz(-2.907739) q[2];
sx q[2];
rz(-0.65143603) q[2];
rz(2.3598119) q[3];
sx q[3];
rz(-1.6580481) q[3];
sx q[3];
rz(-2.206291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5930138) q[0];
sx q[0];
rz(-0.23603708) q[0];
sx q[0];
rz(-0.4655984) q[0];
rz(1.1511401) q[1];
sx q[1];
rz(-1.2460037) q[1];
sx q[1];
rz(1.8021072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27522408) q[0];
sx q[0];
rz(-1.0491519) q[0];
sx q[0];
rz(2.0417968) q[0];
rz(-pi) q[1];
rz(3.1033045) q[2];
sx q[2];
rz(-1.6147505) q[2];
sx q[2];
rz(0.15434855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0995746) q[1];
sx q[1];
rz(-0.23222831) q[1];
sx q[1];
rz(-1.2378511) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0520949) q[3];
sx q[3];
rz(-2.4065402) q[3];
sx q[3];
rz(3.0714761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0307978) q[2];
sx q[2];
rz(-1.75533) q[2];
sx q[2];
rz(-3.1286855) q[2];
rz(-3.0068093) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(0.24136647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4607234) q[0];
sx q[0];
rz(-0.091982059) q[0];
sx q[0];
rz(-1.9665834) q[0];
rz(2.3911632) q[1];
sx q[1];
rz(-1.3576327) q[1];
sx q[1];
rz(-0.39355412) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8581482) q[0];
sx q[0];
rz(-0.99020105) q[0];
sx q[0];
rz(0.20943187) q[0];
rz(-1.4541164) q[2];
sx q[2];
rz(-0.91800729) q[2];
sx q[2];
rz(-2.0443001) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2149512) q[1];
sx q[1];
rz(-1.9991753) q[1];
sx q[1];
rz(0.12081103) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2116488) q[3];
sx q[3];
rz(-2.5287147) q[3];
sx q[3];
rz(1.5604492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5081818) q[2];
sx q[2];
rz(-1.1339302) q[2];
sx q[2];
rz(-0.22937648) q[2];
rz(-0.091382787) q[3];
sx q[3];
rz(-1.4437851) q[3];
sx q[3];
rz(2.5858333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9230187) q[0];
sx q[0];
rz(-2.3526683) q[0];
sx q[0];
rz(-0.58513418) q[0];
rz(1.2773889) q[1];
sx q[1];
rz(-1.2150512) q[1];
sx q[1];
rz(2.9871984) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6638724) q[0];
sx q[0];
rz(-1.1363875) q[0];
sx q[0];
rz(1.4503195) q[0];
rz(-pi) q[1];
rz(-0.67382185) q[2];
sx q[2];
rz(-0.5267064) q[2];
sx q[2];
rz(-0.4835085) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8789492) q[1];
sx q[1];
rz(-0.41987881) q[1];
sx q[1];
rz(2.4559569) q[1];
rz(-pi) q[2];
rz(-1.2326689) q[3];
sx q[3];
rz(-2.6498389) q[3];
sx q[3];
rz(-1.7724747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.020891) q[2];
sx q[2];
rz(-2.2521844) q[2];
sx q[2];
rz(2.6013539) q[2];
rz(0.39515105) q[3];
sx q[3];
rz(-2.4873147) q[3];
sx q[3];
rz(1.3490217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2111135) q[0];
sx q[0];
rz(-1.8607288) q[0];
sx q[0];
rz(-1.6092009) q[0];
rz(0.92597517) q[1];
sx q[1];
rz(-1.4264868) q[1];
sx q[1];
rz(1.6645974) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2944082) q[0];
sx q[0];
rz(-1.5715181) q[0];
sx q[0];
rz(-0.72615726) q[0];
rz(-pi) q[1];
rz(-1.1941998) q[2];
sx q[2];
rz(-1.8459326) q[2];
sx q[2];
rz(-1.9614432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3850579) q[1];
sx q[1];
rz(-1.9137662) q[1];
sx q[1];
rz(-2.1567731) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1891045) q[3];
sx q[3];
rz(-2.455061) q[3];
sx q[3];
rz(-1.8812979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8605139) q[2];
sx q[2];
rz(-1.3382567) q[2];
sx q[2];
rz(0.64858428) q[2];
rz(0.72566882) q[3];
sx q[3];
rz(-0.23894335) q[3];
sx q[3];
rz(1.3919977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7657179) q[0];
sx q[0];
rz(-2.1633042) q[0];
sx q[0];
rz(-2.3194585) q[0];
rz(-2.0769465) q[1];
sx q[1];
rz(-1.6864265) q[1];
sx q[1];
rz(2.0815157) q[1];
rz(2.5049985) q[2];
sx q[2];
rz(-1.7117715) q[2];
sx q[2];
rz(-2.911917) q[2];
rz(-0.81448803) q[3];
sx q[3];
rz(-1.1753488) q[3];
sx q[3];
rz(-1.4594441) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

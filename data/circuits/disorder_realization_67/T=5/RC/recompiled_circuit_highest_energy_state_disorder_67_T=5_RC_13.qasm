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
rz(-1.851097) q[1];
sx q[1];
rz(-1.8654856) q[1];
sx q[1];
rz(-2.2955503) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94933287) q[0];
sx q[0];
rz(-2.7547382) q[0];
sx q[0];
rz(-2.7036376) q[0];
rz(-pi) q[1];
rz(2.8522367) q[2];
sx q[2];
rz(-1.2086007) q[2];
sx q[2];
rz(-1.8487428) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5436094) q[1];
sx q[1];
rz(-2.5777237) q[1];
sx q[1];
rz(-0.51582526) q[1];
x q[2];
rz(-0.11183864) q[3];
sx q[3];
rz(-1.5790644) q[3];
sx q[3];
rz(-1.1638851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35752615) q[2];
sx q[2];
rz(-2.1937723) q[2];
sx q[2];
rz(-0.86304647) q[2];
rz(2.0283902) q[3];
sx q[3];
rz(-2.2200255) q[3];
sx q[3];
rz(-2.34288) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6031826) q[0];
sx q[0];
rz(-0.69284678) q[0];
sx q[0];
rz(-0.26042724) q[0];
rz(0.32568359) q[1];
sx q[1];
rz(-0.98418701) q[1];
sx q[1];
rz(0.30776417) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4935319) q[0];
sx q[0];
rz(-2.2770398) q[0];
sx q[0];
rz(0.34428455) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2529875) q[2];
sx q[2];
rz(-2.0205287) q[2];
sx q[2];
rz(-0.98435246) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5883397) q[1];
sx q[1];
rz(-0.25811895) q[1];
sx q[1];
rz(1.2023167) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.006853718) q[3];
sx q[3];
rz(-0.93539507) q[3];
sx q[3];
rz(1.3069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3088358) q[2];
sx q[2];
rz(-0.86897659) q[2];
sx q[2];
rz(0.55574179) q[2];
rz(-2.6293788) q[3];
sx q[3];
rz(-0.54849505) q[3];
sx q[3];
rz(-0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25278768) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(-0.74932253) q[0];
rz(-2.0454171) q[1];
sx q[1];
rz(-1.7671894) q[1];
sx q[1];
rz(-0.38415092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2040981) q[0];
sx q[0];
rz(-0.9956313) q[0];
sx q[0];
rz(1.8857486) q[0];
rz(1.8245424) q[2];
sx q[2];
rz(-0.98461005) q[2];
sx q[2];
rz(-1.7895928) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.13680116) q[1];
sx q[1];
rz(-0.30711949) q[1];
sx q[1];
rz(1.7476487) q[1];
x q[2];
rz(2.8504653) q[3];
sx q[3];
rz(-2.2343221) q[3];
sx q[3];
rz(-0.71602589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38203794) q[2];
sx q[2];
rz(-1.3052772) q[2];
sx q[2];
rz(1.5142415) q[2];
rz(-2.4115327) q[3];
sx q[3];
rz(-2.4000013) q[3];
sx q[3];
rz(2.7470284) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6877614) q[0];
sx q[0];
rz(-1.7173959) q[0];
sx q[0];
rz(-3.0815531) q[0];
rz(2.2278348) q[1];
sx q[1];
rz(-0.91122952) q[1];
sx q[1];
rz(-2.2027016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37591463) q[0];
sx q[0];
rz(-1.6041099) q[0];
sx q[0];
rz(-0.062996431) q[0];
rz(-pi) q[1];
rz(0.36146116) q[2];
sx q[2];
rz(-0.088080125) q[2];
sx q[2];
rz(0.19780209) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9524433) q[1];
sx q[1];
rz(-1.5464142) q[1];
sx q[1];
rz(-1.0827176) q[1];
x q[2];
rz(-2.3653829) q[3];
sx q[3];
rz(-0.84899711) q[3];
sx q[3];
rz(1.9263903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.011220304) q[2];
sx q[2];
rz(-0.29272407) q[2];
sx q[2];
rz(0.65227738) q[2];
rz(-1.8782714) q[3];
sx q[3];
rz(-1.5247366) q[3];
sx q[3];
rz(1.1711082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0703099) q[0];
sx q[0];
rz(-1.2391397) q[0];
sx q[0];
rz(-2.5791445) q[0];
rz(-1.7930454) q[1];
sx q[1];
rz(-0.28678539) q[1];
sx q[1];
rz(0.58132344) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4084449) q[0];
sx q[0];
rz(-2.3075891) q[0];
sx q[0];
rz(2.187378) q[0];
x q[1];
rz(-2.6810886) q[2];
sx q[2];
rz(-1.5384314) q[2];
sx q[2];
rz(-0.89214395) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2723734) q[1];
sx q[1];
rz(-1.0989184) q[1];
sx q[1];
rz(-2.7414338) q[1];
rz(-pi) q[2];
rz(1.6147037) q[3];
sx q[3];
rz(-0.9287408) q[3];
sx q[3];
rz(1.8521483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.44744667) q[2];
sx q[2];
rz(-1.7116825) q[2];
sx q[2];
rz(-2.4549129) q[2];
rz(0.40003362) q[3];
sx q[3];
rz(-1.8828078) q[3];
sx q[3];
rz(0.011628477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7069063) q[0];
sx q[0];
rz(-0.97395581) q[0];
sx q[0];
rz(0.45027012) q[0];
rz(-0.88092583) q[1];
sx q[1];
rz(-1.2342781) q[1];
sx q[1];
rz(-1.1164104) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4661144) q[0];
sx q[0];
rz(-0.7921392) q[0];
sx q[0];
rz(1.303543) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35043244) q[2];
sx q[2];
rz(-1.3964126) q[2];
sx q[2];
rz(0.36878219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.19710572) q[1];
sx q[1];
rz(-0.61054269) q[1];
sx q[1];
rz(1.6050299) q[1];
x q[2];
rz(-2.0209337) q[3];
sx q[3];
rz(-1.2613457) q[3];
sx q[3];
rz(-1.3510973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66699666) q[2];
sx q[2];
rz(-2.907739) q[2];
sx q[2];
rz(2.4901566) q[2];
rz(-2.3598119) q[3];
sx q[3];
rz(-1.4835446) q[3];
sx q[3];
rz(0.93530161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5485789) q[0];
sx q[0];
rz(-2.9055556) q[0];
sx q[0];
rz(0.4655984) q[0];
rz(-1.1511401) q[1];
sx q[1];
rz(-1.8955889) q[1];
sx q[1];
rz(1.8021072) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5440774) q[0];
sx q[0];
rz(-1.9751514) q[0];
sx q[0];
rz(-0.57283516) q[0];
x q[1];
rz(0.8546245) q[2];
sx q[2];
rz(-0.058283866) q[2];
sx q[2];
rz(2.5789277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20414824) q[1];
sx q[1];
rz(-1.6460858) q[1];
sx q[1];
rz(1.7906899) q[1];
rz(-pi) q[2];
rz(-2.2364113) q[3];
sx q[3];
rz(-1.2318805) q[3];
sx q[3];
rz(1.9012332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0307978) q[2];
sx q[2];
rz(-1.75533) q[2];
sx q[2];
rz(3.1286855) q[2];
rz(0.13478336) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(0.24136647) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4607234) q[0];
sx q[0];
rz(-0.091982059) q[0];
sx q[0];
rz(-1.1750093) q[0];
rz(0.75042945) q[1];
sx q[1];
rz(-1.3576327) q[1];
sx q[1];
rz(0.39355412) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9703044) q[0];
sx q[0];
rz(-1.7455202) q[0];
sx q[0];
rz(2.1615684) q[0];
rz(-2.9904994) q[2];
sx q[2];
rz(-2.4799621) q[2];
sx q[2];
rz(-2.2349295) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4992884) q[1];
sx q[1];
rz(-0.44407108) q[1];
sx q[1];
rz(1.3127692) q[1];
x q[2];
rz(-0.39799798) q[3];
sx q[3];
rz(-2.0500215) q[3];
sx q[3];
rz(-0.82113876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.63341081) q[2];
sx q[2];
rz(-1.1339302) q[2];
sx q[2];
rz(-2.9122162) q[2];
rz(-3.0502099) q[3];
sx q[3];
rz(-1.4437851) q[3];
sx q[3];
rz(-2.5858333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9230187) q[0];
sx q[0];
rz(-2.3526683) q[0];
sx q[0];
rz(0.58513418) q[0];
rz(-1.8642037) q[1];
sx q[1];
rz(-1.9265415) q[1];
sx q[1];
rz(-2.9871984) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0421678) q[0];
sx q[0];
rz(-1.6800361) q[0];
sx q[0];
rz(2.7043998) q[0];
x q[1];
rz(-0.67382185) q[2];
sx q[2];
rz(-0.5267064) q[2];
sx q[2];
rz(-0.4835085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2626434) q[1];
sx q[1];
rz(-2.7217138) q[1];
sx q[1];
rz(-2.4559569) q[1];
x q[2];
rz(2.0386857) q[3];
sx q[3];
rz(-1.7280735) q[3];
sx q[3];
rz(-3.0426971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1207017) q[2];
sx q[2];
rz(-0.88940826) q[2];
sx q[2];
rz(-0.54023877) q[2];
rz(0.39515105) q[3];
sx q[3];
rz(-2.4873147) q[3];
sx q[3];
rz(1.3490217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9304792) q[0];
sx q[0];
rz(-1.8607288) q[0];
sx q[0];
rz(1.6092009) q[0];
rz(2.2156175) q[1];
sx q[1];
rz(-1.4264868) q[1];
sx q[1];
rz(-1.6645974) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4173399) q[0];
sx q[0];
rz(-0.84463929) q[0];
sx q[0];
rz(1.5698311) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2254752) q[2];
sx q[2];
rz(-0.46248702) q[2];
sx q[2];
rz(-0.99257591) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.033869318) q[1];
sx q[1];
rz(-1.0230015) q[1];
sx q[1];
rz(-0.40489991) q[1];
rz(2.6981101) q[3];
sx q[3];
rz(-1.0280307) q[3];
sx q[3];
rz(0.51668985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28107873) q[2];
sx q[2];
rz(-1.8033359) q[2];
sx q[2];
rz(0.64858428) q[2];
rz(-2.4159238) q[3];
sx q[3];
rz(-2.9026493) q[3];
sx q[3];
rz(-1.3919977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3758748) q[0];
sx q[0];
rz(-2.1633042) q[0];
sx q[0];
rz(-2.3194585) q[0];
rz(-2.0769465) q[1];
sx q[1];
rz(-1.6864265) q[1];
sx q[1];
rz(2.0815157) q[1];
rz(-2.5049985) q[2];
sx q[2];
rz(-1.4298212) q[2];
sx q[2];
rz(0.22967568) q[2];
rz(-0.52100427) q[3];
sx q[3];
rz(-0.88501265) q[3];
sx q[3];
rz(2.9043502) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

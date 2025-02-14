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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48117366) q[0];
sx q[0];
rz(-1.9194897) q[0];
sx q[0];
rz(-1.3997188) q[0];
rz(-pi) q[1];
rz(-1.1943075) q[2];
sx q[2];
rz(-1.8408911) q[2];
sx q[2];
rz(0.17284753) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5436094) q[1];
sx q[1];
rz(-0.56386891) q[1];
sx q[1];
rz(-0.51582526) q[1];
x q[2];
rz(-0.073949532) q[3];
sx q[3];
rz(-3.0294501) q[3];
sx q[3];
rz(-2.6611947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7840665) q[2];
sx q[2];
rz(-2.1937723) q[2];
sx q[2];
rz(0.86304647) q[2];
rz(2.0283902) q[3];
sx q[3];
rz(-0.92156711) q[3];
sx q[3];
rz(-0.79871261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5384101) q[0];
sx q[0];
rz(-2.4487459) q[0];
sx q[0];
rz(2.8811654) q[0];
rz(2.8159091) q[1];
sx q[1];
rz(-2.1574056) q[1];
sx q[1];
rz(0.30776417) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4474898) q[0];
sx q[0];
rz(-1.3110975) q[0];
sx q[0];
rz(-0.83456852) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5670916) q[2];
sx q[2];
rz(-0.5443474) q[2];
sx q[2];
rz(-0.33657227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.766749) q[1];
sx q[1];
rz(-1.6628712) q[1];
sx q[1];
rz(1.8122775) q[1];
x q[2];
rz(-3.1347389) q[3];
sx q[3];
rz(-0.93539507) q[3];
sx q[3];
rz(1.8346661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3088358) q[2];
sx q[2];
rz(-2.2726161) q[2];
sx q[2];
rz(0.55574179) q[2];
rz(2.6293788) q[3];
sx q[3];
rz(-2.5930976) q[3];
sx q[3];
rz(-0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.888805) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(2.3922701) q[0];
rz(2.0454171) q[1];
sx q[1];
rz(-1.3744033) q[1];
sx q[1];
rz(2.7574417) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644728) q[0];
sx q[0];
rz(-0.64711231) q[0];
sx q[0];
rz(-0.44575925) q[0];
x q[1];
rz(-0.60127705) q[2];
sx q[2];
rz(-1.3601175) q[2];
sx q[2];
rz(-0.36128584) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2652581) q[1];
sx q[1];
rz(-1.6240082) q[1];
sx q[1];
rz(1.8734147) q[1];
rz(2.2552842) q[3];
sx q[3];
rz(-1.7989007) q[3];
sx q[3];
rz(-2.4693173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7595547) q[2];
sx q[2];
rz(-1.3052772) q[2];
sx q[2];
rz(1.6273512) q[2];
rz(0.73005992) q[3];
sx q[3];
rz(-0.74159139) q[3];
sx q[3];
rz(0.39456427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.45383129) q[0];
sx q[0];
rz(-1.4241968) q[0];
sx q[0];
rz(0.060039595) q[0];
rz(-0.91375786) q[1];
sx q[1];
rz(-2.2303631) q[1];
sx q[1];
rz(-0.93889108) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9446099) q[0];
sx q[0];
rz(-1.6337578) q[0];
sx q[0];
rz(1.5374166) q[0];
rz(-pi) q[1];
rz(-2.7801315) q[2];
sx q[2];
rz(-0.088080125) q[2];
sx q[2];
rz(-2.9437906) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7728887) q[1];
sx q[1];
rz(-2.0587173) q[1];
sx q[1];
rz(-0.027603961) q[1];
rz(-0.89857159) q[3];
sx q[3];
rz(-1.0055526) q[3];
sx q[3];
rz(2.1936301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1303723) q[2];
sx q[2];
rz(-2.8488686) q[2];
sx q[2];
rz(-0.65227738) q[2];
rz(1.8782714) q[3];
sx q[3];
rz(-1.5247366) q[3];
sx q[3];
rz(-1.1711082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0703099) q[0];
sx q[0];
rz(-1.2391397) q[0];
sx q[0];
rz(2.5791445) q[0];
rz(1.7930454) q[1];
sx q[1];
rz(-2.8548073) q[1];
sx q[1];
rz(0.58132344) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5347915) q[0];
sx q[0];
rz(-1.1282217) q[0];
sx q[0];
rz(2.3032196) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.072725459) q[2];
sx q[2];
rz(-0.46155864) q[2];
sx q[2];
rz(2.3978021) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6186468) q[1];
sx q[1];
rz(-0.60877555) q[1];
sx q[1];
rz(2.2227915) q[1];
rz(1.526889) q[3];
sx q[3];
rz(-0.9287408) q[3];
sx q[3];
rz(1.2894443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44744667) q[2];
sx q[2];
rz(-1.7116825) q[2];
sx q[2];
rz(0.68667975) q[2];
rz(2.741559) q[3];
sx q[3];
rz(-1.8828078) q[3];
sx q[3];
rz(3.1299642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7069063) q[0];
sx q[0];
rz(-0.97395581) q[0];
sx q[0];
rz(0.45027012) q[0];
rz(0.88092583) q[1];
sx q[1];
rz(-1.2342781) q[1];
sx q[1];
rz(-2.0251822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6754782) q[0];
sx q[0];
rz(-0.7921392) q[0];
sx q[0];
rz(1.8380497) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7562148) q[2];
sx q[2];
rz(-1.2259019) q[2];
sx q[2];
rz(-1.138681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3456381) q[1];
sx q[1];
rz(-1.5511723) q[1];
sx q[1];
rz(0.96052891) q[1];
x q[2];
rz(-2.2045256) q[3];
sx q[3];
rz(-0.54017276) q[3];
sx q[3];
rz(-0.34264178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66699666) q[2];
sx q[2];
rz(-0.23385364) q[2];
sx q[2];
rz(-2.4901566) q[2];
rz(0.78178072) q[3];
sx q[3];
rz(-1.6580481) q[3];
sx q[3];
rz(2.206291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5485789) q[0];
sx q[0];
rz(-0.23603708) q[0];
sx q[0];
rz(-2.6759942) q[0];
rz(1.1511401) q[1];
sx q[1];
rz(-1.2460037) q[1];
sx q[1];
rz(-1.3394855) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5975153) q[0];
sx q[0];
rz(-1.1664412) q[0];
sx q[0];
rz(0.57283516) q[0];
rz(3.1033045) q[2];
sx q[2];
rz(-1.5268421) q[2];
sx q[2];
rz(-0.15434855) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9374444) q[1];
sx q[1];
rz(-1.4955069) q[1];
sx q[1];
rz(1.3509028) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0894977) q[3];
sx q[3];
rz(-0.7350525) q[3];
sx q[3];
rz(0.07011655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1107948) q[2];
sx q[2];
rz(-1.3862627) q[2];
sx q[2];
rz(-0.01290713) q[2];
rz(3.0068093) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(2.9002262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1712883) q[0];
sx q[0];
rz(-1.3960724) q[0];
sx q[0];
rz(0.98002429) q[0];
rz(-pi) q[1];
rz(-0.15109328) q[2];
sx q[2];
rz(-2.4799621) q[2];
sx q[2];
rz(-0.90666319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8361349) q[1];
sx q[1];
rz(-1.6806446) q[1];
sx q[1];
rz(-1.1396465) q[1];
rz(-pi) q[2];
rz(-2.7435947) q[3];
sx q[3];
rz(-2.0500215) q[3];
sx q[3];
rz(-2.3204539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5081818) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9230187) q[0];
sx q[0];
rz(-0.78892437) q[0];
sx q[0];
rz(-2.5564585) q[0];
rz(1.8642037) q[1];
sx q[1];
rz(-1.2150512) q[1];
sx q[1];
rz(0.1543943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6638724) q[0];
sx q[0];
rz(-2.0052052) q[0];
sx q[0];
rz(-1.4503195) q[0];
rz(-pi) q[1];
rz(0.67382185) q[2];
sx q[2];
rz(-0.5267064) q[2];
sx q[2];
rz(-2.6580842) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8789492) q[1];
sx q[1];
rz(-2.7217138) q[1];
sx q[1];
rz(-2.4559569) q[1];
rz(-pi) q[2];
rz(-2.9657439) q[3];
sx q[3];
rz(-1.1091341) q[3];
sx q[3];
rz(-1.3929183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.020891) q[2];
sx q[2];
rz(-2.2521844) q[2];
sx q[2];
rz(-2.6013539) q[2];
rz(-2.7464416) q[3];
sx q[3];
rz(-2.4873147) q[3];
sx q[3];
rz(1.3490217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2111135) q[0];
sx q[0];
rz(-1.8607288) q[0];
sx q[0];
rz(1.6092009) q[0];
rz(2.2156175) q[1];
sx q[1];
rz(-1.7151058) q[1];
sx q[1];
rz(-1.4769953) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2944082) q[0];
sx q[0];
rz(-1.5715181) q[0];
sx q[0];
rz(0.72615726) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2254752) q[2];
sx q[2];
rz(-2.6791056) q[2];
sx q[2];
rz(-0.99257591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3850579) q[1];
sx q[1];
rz(-1.9137662) q[1];
sx q[1];
rz(-0.98481957) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44348251) q[3];
sx q[3];
rz(-1.0280307) q[3];
sx q[3];
rz(-0.51668985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28107873) q[2];
sx q[2];
rz(-1.3382567) q[2];
sx q[2];
rz(2.4930084) q[2];
rz(-0.72566882) q[3];
sx q[3];
rz(-2.9026493) q[3];
sx q[3];
rz(1.3919977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7657179) q[0];
sx q[0];
rz(-2.1633042) q[0];
sx q[0];
rz(-2.3194585) q[0];
rz(2.0769465) q[1];
sx q[1];
rz(-1.4551661) q[1];
sx q[1];
rz(-1.0600769) q[1];
rz(-1.3961095) q[2];
sx q[2];
rz(-2.2000763) q[2];
sx q[2];
rz(1.6969703) q[2];
rz(2.1172932) q[3];
sx q[3];
rz(-0.83491769) q[3];
sx q[3];
rz(0.49900311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

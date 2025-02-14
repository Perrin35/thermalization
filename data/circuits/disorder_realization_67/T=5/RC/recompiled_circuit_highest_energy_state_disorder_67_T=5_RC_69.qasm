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
rz(1.2904957) q[1];
sx q[1];
rz(-1.2761071) q[1];
sx q[1];
rz(-0.84604231) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94933287) q[0];
sx q[0];
rz(-2.7547382) q[0];
sx q[0];
rz(0.43795507) q[0];
x q[1];
rz(-2.2162391) q[2];
sx q[2];
rz(-2.6820289) q[2];
sx q[2];
rz(-1.1499576) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5979833) q[1];
sx q[1];
rz(-0.56386891) q[1];
sx q[1];
rz(-2.6257674) q[1];
x q[2];
rz(0.11183864) q[3];
sx q[3];
rz(-1.5625283) q[3];
sx q[3];
rz(-1.1638851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7840665) q[2];
sx q[2];
rz(-0.94782031) q[2];
sx q[2];
rz(0.86304647) q[2];
rz(-1.1132025) q[3];
sx q[3];
rz(-2.2200255) q[3];
sx q[3];
rz(0.79871261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5384101) q[0];
sx q[0];
rz(-2.4487459) q[0];
sx q[0];
rz(-0.26042724) q[0];
rz(-2.8159091) q[1];
sx q[1];
rz(-0.98418701) q[1];
sx q[1];
rz(0.30776417) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1528252) q[0];
sx q[0];
rz(-0.77251311) q[0];
sx q[0];
rz(1.9475742) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5670916) q[2];
sx q[2];
rz(-2.5972453) q[2];
sx q[2];
rz(-0.33657227) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.766749) q[1];
sx q[1];
rz(-1.6628712) q[1];
sx q[1];
rz(-1.3293152) q[1];
rz(-0.93538385) q[3];
sx q[3];
rz(-1.5652802) q[3];
sx q[3];
rz(2.8817906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8327568) q[2];
sx q[2];
rz(-2.2726161) q[2];
sx q[2];
rz(0.55574179) q[2];
rz(-2.6293788) q[3];
sx q[3];
rz(-2.5930976) q[3];
sx q[3];
rz(0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.888805) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(2.3922701) q[0];
rz(-1.0961756) q[1];
sx q[1];
rz(-1.3744033) q[1];
sx q[1];
rz(2.7574417) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4771198) q[0];
sx q[0];
rz(-2.4944803) q[0];
sx q[0];
rz(-2.6958334) q[0];
rz(1.3170502) q[2];
sx q[2];
rz(-2.1569826) q[2];
sx q[2];
rz(-1.7895928) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0047915) q[1];
sx q[1];
rz(-0.30711949) q[1];
sx q[1];
rz(-1.3939439) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88630845) q[3];
sx q[3];
rz(-1.3426919) q[3];
sx q[3];
rz(0.67227539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38203794) q[2];
sx q[2];
rz(-1.8363154) q[2];
sx q[2];
rz(-1.5142415) q[2];
rz(-2.4115327) q[3];
sx q[3];
rz(-0.74159139) q[3];
sx q[3];
rz(0.39456427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6877614) q[0];
sx q[0];
rz(-1.4241968) q[0];
sx q[0];
rz(3.0815531) q[0];
rz(-0.91375786) q[1];
sx q[1];
rz(-0.91122952) q[1];
sx q[1];
rz(0.93889108) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37591463) q[0];
sx q[0];
rz(-1.5374827) q[0];
sx q[0];
rz(3.0785962) q[0];
rz(2.7801315) q[2];
sx q[2];
rz(-0.088080125) q[2];
sx q[2];
rz(2.9437906) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42753427) q[1];
sx q[1];
rz(-0.48863829) q[1];
sx q[1];
rz(-1.6227551) q[1];
rz(-pi) q[2];
rz(2.3653829) q[3];
sx q[3];
rz(-2.2925955) q[3];
sx q[3];
rz(-1.2152023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1303723) q[2];
sx q[2];
rz(-0.29272407) q[2];
sx q[2];
rz(-2.4893153) q[2];
rz(-1.8782714) q[3];
sx q[3];
rz(-1.616856) q[3];
sx q[3];
rz(-1.1711082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071282722) q[0];
sx q[0];
rz(-1.902453) q[0];
sx q[0];
rz(-0.5624482) q[0];
rz(1.7930454) q[1];
sx q[1];
rz(-0.28678539) q[1];
sx q[1];
rz(2.5602692) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59636694) q[0];
sx q[0];
rz(-2.2196182) q[0];
sx q[0];
rz(-0.56745402) q[0];
x q[1];
rz(-0.072725459) q[2];
sx q[2];
rz(-2.680034) q[2];
sx q[2];
rz(0.74379057) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86921924) q[1];
sx q[1];
rz(-2.0426742) q[1];
sx q[1];
rz(-2.7414338) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6147037) q[3];
sx q[3];
rz(-2.2128519) q[3];
sx q[3];
rz(1.2894443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.7069063) q[0];
sx q[0];
rz(-2.1676368) q[0];
sx q[0];
rz(2.6913225) q[0];
rz(-2.2606668) q[1];
sx q[1];
rz(-1.2342781) q[1];
sx q[1];
rz(-2.0251822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4661144) q[0];
sx q[0];
rz(-2.3494534) q[0];
sx q[0];
rz(1.8380497) q[0];
rz(-0.35043244) q[2];
sx q[2];
rz(-1.3964126) q[2];
sx q[2];
rz(-2.7728105) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7959545) q[1];
sx q[1];
rz(-1.5511723) q[1];
sx q[1];
rz(-0.96052891) q[1];
x q[2];
rz(-2.2045256) q[3];
sx q[3];
rz(-0.54017276) q[3];
sx q[3];
rz(-0.34264178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66699666) q[2];
sx q[2];
rz(-0.23385364) q[2];
sx q[2];
rz(0.65143603) q[2];
rz(-2.3598119) q[3];
sx q[3];
rz(-1.4835446) q[3];
sx q[3];
rz(0.93530161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5930138) q[0];
sx q[0];
rz(-0.23603708) q[0];
sx q[0];
rz(2.6759942) q[0];
rz(1.1511401) q[1];
sx q[1];
rz(-1.8955889) q[1];
sx q[1];
rz(1.3394855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6205821) q[0];
sx q[0];
rz(-0.68787439) q[0];
sx q[0];
rz(-2.4732711) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2869682) q[2];
sx q[2];
rz(-0.058283866) q[2];
sx q[2];
rz(2.5789277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9374444) q[1];
sx q[1];
rz(-1.6460858) q[1];
sx q[1];
rz(1.3509028) q[1];
rz(-pi) q[2];
rz(-2.0894977) q[3];
sx q[3];
rz(-2.4065402) q[3];
sx q[3];
rz(-3.0714761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0307978) q[2];
sx q[2];
rz(-1.3862627) q[2];
sx q[2];
rz(-0.01290713) q[2];
rz(3.0068093) q[3];
sx q[3];
rz(-1.8811992) q[3];
sx q[3];
rz(-2.9002262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4607234) q[0];
sx q[0];
rz(-3.0496106) q[0];
sx q[0];
rz(1.9665834) q[0];
rz(-0.75042945) q[1];
sx q[1];
rz(-1.3576327) q[1];
sx q[1];
rz(-0.39355412) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9703044) q[0];
sx q[0];
rz(-1.3960724) q[0];
sx q[0];
rz(-0.98002429) q[0];
x q[1];
rz(2.9904994) q[2];
sx q[2];
rz(-2.4799621) q[2];
sx q[2];
rz(-0.90666319) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30545772) q[1];
sx q[1];
rz(-1.4609481) q[1];
sx q[1];
rz(2.0019462) q[1];
rz(-pi) q[2];
rz(2.2116488) q[3];
sx q[3];
rz(-2.5287147) q[3];
sx q[3];
rz(1.5811435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.9230187) q[0];
sx q[0];
rz(-2.3526683) q[0];
sx q[0];
rz(-0.58513418) q[0];
rz(-1.2773889) q[1];
sx q[1];
rz(-1.9265415) q[1];
sx q[1];
rz(-0.1543943) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0421678) q[0];
sx q[0];
rz(-1.4615566) q[0];
sx q[0];
rz(-2.7043998) q[0];
x q[1];
rz(-1.9188668) q[2];
sx q[2];
rz(-1.1670975) q[2];
sx q[2];
rz(-0.26223809) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8789492) q[1];
sx q[1];
rz(-2.7217138) q[1];
sx q[1];
rz(-2.4559569) q[1];
rz(-pi) q[2];
rz(2.9657439) q[3];
sx q[3];
rz(-2.0324586) q[3];
sx q[3];
rz(1.7486743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1207017) q[2];
sx q[2];
rz(-0.88940826) q[2];
sx q[2];
rz(-0.54023877) q[2];
rz(-2.7464416) q[3];
sx q[3];
rz(-2.4873147) q[3];
sx q[3];
rz(-1.7925709) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2111135) q[0];
sx q[0];
rz(-1.8607288) q[0];
sx q[0];
rz(-1.6092009) q[0];
rz(-0.92597517) q[1];
sx q[1];
rz(-1.4264868) q[1];
sx q[1];
rz(-1.6645974) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72425276) q[0];
sx q[0];
rz(-2.2969534) q[0];
sx q[0];
rz(-1.5717616) q[0];
x q[1];
rz(1.9473929) q[2];
sx q[2];
rz(-1.8459326) q[2];
sx q[2];
rz(-1.9614432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1077233) q[1];
sx q[1];
rz(-2.1185912) q[1];
sx q[1];
rz(0.40489991) q[1];
x q[2];
rz(2.1891045) q[3];
sx q[3];
rz(-2.455061) q[3];
sx q[3];
rz(1.8812979) q[3];
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
rz(2.4159238) q[3];
sx q[3];
rz(-0.23894335) q[3];
sx q[3];
rz(1.749595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9072472) q[2];
sx q[2];
rz(-0.6498944) q[2];
sx q[2];
rz(1.9882974) q[2];
rz(-1.0242994) q[3];
sx q[3];
rz(-0.83491769) q[3];
sx q[3];
rz(0.49900311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

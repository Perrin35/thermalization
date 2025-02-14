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
rz(0.30581185) q[0];
sx q[0];
rz(3.7753063) q[0];
sx q[0];
rz(10.033327) q[0];
rz(-0.67148709) q[1];
sx q[1];
rz(-2.4784213) q[1];
sx q[1];
rz(1.4594249) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7072351) q[0];
sx q[0];
rz(-0.57091129) q[0];
sx q[0];
rz(2.0152604) q[0];
rz(2.9179108) q[2];
sx q[2];
rz(-2.0111901) q[2];
sx q[2];
rz(0.5612095) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83460036) q[1];
sx q[1];
rz(-1.9687593) q[1];
sx q[1];
rz(2.2521521) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18078928) q[3];
sx q[3];
rz(-0.76813625) q[3];
sx q[3];
rz(0.013313499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3964316) q[2];
sx q[2];
rz(-0.67239422) q[2];
sx q[2];
rz(-2.3154837) q[2];
rz(2.7769026) q[3];
sx q[3];
rz(-1.5267905) q[3];
sx q[3];
rz(-0.29432347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.9593338) q[0];
sx q[0];
rz(-1.2428281) q[0];
sx q[0];
rz(0.44573319) q[0];
rz(-0.590473) q[1];
sx q[1];
rz(-2.1069374) q[1];
sx q[1];
rz(-2.7139434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2954461) q[0];
sx q[0];
rz(-1.5238191) q[0];
sx q[0];
rz(-2.0278005) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81645963) q[2];
sx q[2];
rz(-2.0000946) q[2];
sx q[2];
rz(-0.72390899) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1688437) q[1];
sx q[1];
rz(-2.4764937) q[1];
sx q[1];
rz(-1.14023) q[1];
rz(-0.79645653) q[3];
sx q[3];
rz(-1.5672308) q[3];
sx q[3];
rz(-2.2623747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4599956) q[2];
sx q[2];
rz(-0.57969105) q[2];
sx q[2];
rz(1.0350636) q[2];
rz(1.5929818) q[3];
sx q[3];
rz(-2.9620453) q[3];
sx q[3];
rz(-2.2523994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0276133) q[0];
sx q[0];
rz(-0.005682156) q[0];
sx q[0];
rz(-3.0300544) q[0];
rz(-1.5533718) q[1];
sx q[1];
rz(-2.7370743) q[1];
sx q[1];
rz(2.9389971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3521393) q[0];
sx q[0];
rz(-1.4367665) q[0];
sx q[0];
rz(1.9699957) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8526788) q[2];
sx q[2];
rz(-1.2855069) q[2];
sx q[2];
rz(2.4288175) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93567551) q[1];
sx q[1];
rz(-0.96921235) q[1];
sx q[1];
rz(2.4755073) q[1];
rz(-pi) q[2];
rz(1.1502017) q[3];
sx q[3];
rz(-2.1199825) q[3];
sx q[3];
rz(-2.7830916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4488039) q[2];
sx q[2];
rz(-2.0710129) q[2];
sx q[2];
rz(0.98133522) q[2];
rz(2.8547309) q[3];
sx q[3];
rz(-2.562264) q[3];
sx q[3];
rz(-2.9441492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9051179) q[0];
sx q[0];
rz(-3.0906257) q[0];
sx q[0];
rz(0.68848759) q[0];
rz(2.9972637) q[1];
sx q[1];
rz(-1.1878443) q[1];
sx q[1];
rz(-0.7483288) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8546974) q[0];
sx q[0];
rz(-2.4208491) q[0];
sx q[0];
rz(3.0740644) q[0];
rz(0.26864895) q[2];
sx q[2];
rz(-1.2792818) q[2];
sx q[2];
rz(2.6644601) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6627745) q[1];
sx q[1];
rz(-1.2726901) q[1];
sx q[1];
rz(-1.7653905) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9311936) q[3];
sx q[3];
rz(-0.73355094) q[3];
sx q[3];
rz(-1.4195132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0904842) q[2];
sx q[2];
rz(-1.8277617) q[2];
sx q[2];
rz(-1.2468106) q[2];
rz(3.0221853) q[3];
sx q[3];
rz(-1.2179008) q[3];
sx q[3];
rz(-2.8363805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3731821) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(-0.61491948) q[0];
rz(-0.83456314) q[1];
sx q[1];
rz(-1.395023) q[1];
sx q[1];
rz(0.83647299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1008719) q[0];
sx q[0];
rz(-0.03558579) q[0];
sx q[0];
rz(1.3077223) q[0];
rz(-0.66906382) q[2];
sx q[2];
rz(-2.0537801) q[2];
sx q[2];
rz(1.3945182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53764105) q[1];
sx q[1];
rz(-1.0102235) q[1];
sx q[1];
rz(-0.13092299) q[1];
rz(-2.3990666) q[3];
sx q[3];
rz(-0.98127413) q[3];
sx q[3];
rz(1.3878915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5323083) q[2];
sx q[2];
rz(-0.60695761) q[2];
sx q[2];
rz(0.92030805) q[2];
rz(-2.3406384) q[3];
sx q[3];
rz(-2.6194173) q[3];
sx q[3];
rz(-0.33867684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1717218) q[0];
sx q[0];
rz(-2.0750177) q[0];
sx q[0];
rz(0.26611662) q[0];
rz(1.7099821) q[1];
sx q[1];
rz(-1.1191198) q[1];
sx q[1];
rz(-2.5439579) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1475567) q[0];
sx q[0];
rz(-1.3809675) q[0];
sx q[0];
rz(-2.5231383) q[0];
rz(1.581988) q[2];
sx q[2];
rz(-1.7518242) q[2];
sx q[2];
rz(-2.343296) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3417334) q[1];
sx q[1];
rz(-1.6266154) q[1];
sx q[1];
rz(2.4444405) q[1];
rz(-pi) q[2];
rz(-1.131874) q[3];
sx q[3];
rz(-2.7898129) q[3];
sx q[3];
rz(1.1949415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1282244) q[2];
sx q[2];
rz(-2.8743663) q[2];
sx q[2];
rz(-0.5989778) q[2];
rz(-2.525575) q[3];
sx q[3];
rz(-0.44121656) q[3];
sx q[3];
rz(1.0429355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1004341) q[0];
sx q[0];
rz(-2.2293595) q[0];
sx q[0];
rz(2.8225733) q[0];
rz(3.0306385) q[1];
sx q[1];
rz(-1.7496505) q[1];
sx q[1];
rz(0.14161938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.148478) q[0];
sx q[0];
rz(-1.6671868) q[0];
sx q[0];
rz(-2.4575811) q[0];
rz(1.5779823) q[2];
sx q[2];
rz(-0.55139667) q[2];
sx q[2];
rz(-1.0137272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31894058) q[1];
sx q[1];
rz(-1.3825065) q[1];
sx q[1];
rz(-0.55231185) q[1];
x q[2];
rz(-1.3124489) q[3];
sx q[3];
rz(-0.53505361) q[3];
sx q[3];
rz(-2.759397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5854599) q[2];
sx q[2];
rz(-1.0254878) q[2];
sx q[2];
rz(0.11036135) q[2];
rz(0.50802556) q[3];
sx q[3];
rz(-0.042162687) q[3];
sx q[3];
rz(-1.149811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0069649) q[0];
sx q[0];
rz(-0.65161324) q[0];
sx q[0];
rz(1.2192669) q[0];
rz(3.1189175) q[1];
sx q[1];
rz(-2.2494648) q[1];
sx q[1];
rz(2.5882744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39529453) q[0];
sx q[0];
rz(-0.91522258) q[0];
sx q[0];
rz(-3.0297999) q[0];
rz(-pi) q[1];
rz(-0.99626096) q[2];
sx q[2];
rz(-1.8991514) q[2];
sx q[2];
rz(2.7789214) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0670315) q[1];
sx q[1];
rz(-2.0843049) q[1];
sx q[1];
rz(-1.1719955) q[1];
rz(2.8944772) q[3];
sx q[3];
rz(-1.860642) q[3];
sx q[3];
rz(0.71699504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7360709) q[2];
sx q[2];
rz(-0.61554474) q[2];
sx q[2];
rz(1.7919398) q[2];
rz(-0.060529709) q[3];
sx q[3];
rz(-2.2318201) q[3];
sx q[3];
rz(2.5133666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26199207) q[0];
sx q[0];
rz(-2.1475726) q[0];
sx q[0];
rz(-3.1349728) q[0];
rz(-2.5026542) q[1];
sx q[1];
rz(-2.9220351) q[1];
sx q[1];
rz(-2.7117597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.854786) q[0];
sx q[0];
rz(-1.9298975) q[0];
sx q[0];
rz(-0.39639985) q[0];
rz(1.1388135) q[2];
sx q[2];
rz(-1.3393992) q[2];
sx q[2];
rz(-2.9803986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66664106) q[1];
sx q[1];
rz(-0.78814298) q[1];
sx q[1];
rz(1.5402543) q[1];
rz(-0.58542975) q[3];
sx q[3];
rz(-1.5133891) q[3];
sx q[3];
rz(1.9229744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2868353) q[2];
sx q[2];
rz(-1.9660051) q[2];
sx q[2];
rz(-0.01469928) q[2];
rz(-2.0059351) q[3];
sx q[3];
rz(-2.0292) q[3];
sx q[3];
rz(1.8358102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056674615) q[0];
sx q[0];
rz(-0.25923964) q[0];
sx q[0];
rz(-1.8490476) q[0];
rz(0.1560642) q[1];
sx q[1];
rz(-2.3239457) q[1];
sx q[1];
rz(-2.3110716) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9915402) q[0];
sx q[0];
rz(-0.73313289) q[0];
sx q[0];
rz(1.7197648) q[0];
rz(-pi) q[1];
rz(-0.79367001) q[2];
sx q[2];
rz(-1.7883232) q[2];
sx q[2];
rz(0.95784369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76354181) q[1];
sx q[1];
rz(-1.6507859) q[1];
sx q[1];
rz(-1.3536118) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51955207) q[3];
sx q[3];
rz(-0.41659912) q[3];
sx q[3];
rz(-2.8263457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.10892756) q[2];
sx q[2];
rz(-2.445161) q[2];
sx q[2];
rz(2.6984974) q[2];
rz(-0.14919925) q[3];
sx q[3];
rz(-0.34330338) q[3];
sx q[3];
rz(0.33299115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14015848) q[0];
sx q[0];
rz(-2.5509111) q[0];
sx q[0];
rz(2.1668707) q[0];
rz(-0.99826605) q[1];
sx q[1];
rz(-1.2405735) q[1];
sx q[1];
rz(2.1539198) q[1];
rz(-1.7723631) q[2];
sx q[2];
rz(-1.2342296) q[2];
sx q[2];
rz(-0.078757503) q[2];
rz(2.2343288) q[3];
sx q[3];
rz(-1.2259007) q[3];
sx q[3];
rz(-2.905734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-2.5295244) q[0];
sx q[0];
rz(3.0906313) q[0];
rz(-1.851097) q[1];
sx q[1];
rz(-1.8654856) q[1];
sx q[1];
rz(-2.2955503) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48117366) q[0];
sx q[0];
rz(-1.2221029) q[0];
sx q[0];
rz(-1.7418738) q[0];
rz(-pi) q[1];
rz(-2.8522367) q[2];
sx q[2];
rz(-1.2086007) q[2];
sx q[2];
rz(1.8487428) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4197246) q[1];
sx q[1];
rz(-1.3040191) q[1];
sx q[1];
rz(2.6386923) q[1];
rz(-pi) q[2];
rz(0.11183864) q[3];
sx q[3];
rz(-1.5790644) q[3];
sx q[3];
rz(1.1638851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35752615) q[2];
sx q[2];
rz(-2.1937723) q[2];
sx q[2];
rz(-2.2785462) q[2];
rz(1.1132025) q[3];
sx q[3];
rz(-2.2200255) q[3];
sx q[3];
rz(2.34288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5384101) q[0];
sx q[0];
rz(-2.4487459) q[0];
sx q[0];
rz(-0.26042724) q[0];
rz(2.8159091) q[1];
sx q[1];
rz(-2.1574056) q[1];
sx q[1];
rz(-2.8338285) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1528252) q[0];
sx q[0];
rz(-0.77251311) q[0];
sx q[0];
rz(1.9475742) q[0];
rz(-2.671428) q[2];
sx q[2];
rz(-1.2855296) q[2];
sx q[2];
rz(0.72848749) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1733117) q[1];
sx q[1];
rz(-1.8112343) q[1];
sx q[1];
rz(0.094810061) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2062088) q[3];
sx q[3];
rz(-1.5652802) q[3];
sx q[3];
rz(-2.8817906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3088358) q[2];
sx q[2];
rz(-0.86897659) q[2];
sx q[2];
rz(-0.55574179) q[2];
rz(-0.51221383) q[3];
sx q[3];
rz(-0.54849505) q[3];
sx q[3];
rz(-2.6223124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25278768) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(2.3922701) q[0];
rz(2.0454171) q[1];
sx q[1];
rz(-1.7671894) q[1];
sx q[1];
rz(-2.7574417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45789832) q[0];
sx q[0];
rz(-1.3078469) q[0];
sx q[0];
rz(-2.5431741) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3170502) q[2];
sx q[2];
rz(-2.1569826) q[2];
sx q[2];
rz(-1.7895928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8763346) q[1];
sx q[1];
rz(-1.5175845) q[1];
sx q[1];
rz(-1.268178) q[1];
rz(-pi) q[2];
rz(-2.2552842) q[3];
sx q[3];
rz(-1.7989007) q[3];
sx q[3];
rz(-0.67227539) q[3];
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
rz(2.4115327) q[3];
sx q[3];
rz(-0.74159139) q[3];
sx q[3];
rz(2.7470284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45383129) q[0];
sx q[0];
rz(-1.4241968) q[0];
sx q[0];
rz(-0.060039595) q[0];
rz(0.91375786) q[1];
sx q[1];
rz(-0.91122952) q[1];
sx q[1];
rz(2.2027016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9446099) q[0];
sx q[0];
rz(-1.6337578) q[0];
sx q[0];
rz(-1.5374166) q[0];
x q[1];
rz(1.6020158) q[2];
sx q[2];
rz(-1.4884212) q[2];
sx q[2];
rz(0.16494575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1891494) q[1];
sx q[1];
rz(-1.5951785) q[1];
sx q[1];
rz(2.0588751) q[1];
x q[2];
rz(-2.3653829) q[3];
sx q[3];
rz(-2.2925955) q[3];
sx q[3];
rz(-1.9263903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1303723) q[2];
sx q[2];
rz(-0.29272407) q[2];
sx q[2];
rz(-2.4893153) q[2];
rz(1.8782714) q[3];
sx q[3];
rz(-1.5247366) q[3];
sx q[3];
rz(-1.1711082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071282722) q[0];
sx q[0];
rz(-1.2391397) q[0];
sx q[0];
rz(2.5791445) q[0];
rz(-1.7930454) q[1];
sx q[1];
rz(-0.28678539) q[1];
sx q[1];
rz(0.58132344) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7331478) q[0];
sx q[0];
rz(-0.83400351) q[0];
sx q[0];
rz(-0.95421469) q[0];
rz(0.46050408) q[2];
sx q[2];
rz(-1.5384314) q[2];
sx q[2];
rz(2.2494487) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5229458) q[1];
sx q[1];
rz(-2.5328171) q[1];
sx q[1];
rz(-2.2227915) q[1];
rz(0.64251803) q[3];
sx q[3];
rz(-1.6059562) q[3];
sx q[3];
rz(-2.8865451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.694146) q[2];
sx q[2];
rz(-1.7116825) q[2];
sx q[2];
rz(-2.4549129) q[2];
rz(0.40003362) q[3];
sx q[3];
rz(-1.2587849) q[3];
sx q[3];
rz(-0.011628477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.43468633) q[0];
sx q[0];
rz(-0.97395581) q[0];
sx q[0];
rz(0.45027012) q[0];
rz(-2.2606668) q[1];
sx q[1];
rz(-1.2342781) q[1];
sx q[1];
rz(-2.0251822) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8378431) q[0];
sx q[0];
rz(-0.81401304) q[0];
sx q[0];
rz(-0.26153691) q[0];
rz(2.7911602) q[2];
sx q[2];
rz(-1.74518) q[2];
sx q[2];
rz(-0.36878219) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.23887983) q[1];
sx q[1];
rz(-2.1809291) q[1];
sx q[1];
rz(3.1176477) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1206589) q[3];
sx q[3];
rz(-1.2613457) q[3];
sx q[3];
rz(1.3510973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.474596) q[2];
sx q[2];
rz(-0.23385364) q[2];
sx q[2];
rz(0.65143603) q[2];
rz(2.3598119) q[3];
sx q[3];
rz(-1.4835446) q[3];
sx q[3];
rz(2.206291) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5975153) q[0];
sx q[0];
rz(-1.9751514) q[0];
sx q[0];
rz(-2.5687575) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8546245) q[2];
sx q[2];
rz(-0.058283866) q[2];
sx q[2];
rz(0.56266498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0995746) q[1];
sx q[1];
rz(-2.9093643) q[1];
sx q[1];
rz(-1.2378511) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0520949) q[3];
sx q[3];
rz(-0.7350525) q[3];
sx q[3];
rz(3.0714761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0307978) q[2];
sx q[2];
rz(-1.75533) q[2];
sx q[2];
rz(-0.01290713) q[2];
rz(-3.0068093) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(-2.9002262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4607234) q[0];
sx q[0];
rz(-3.0496106) q[0];
sx q[0];
rz(1.1750093) q[0];
rz(0.75042945) q[1];
sx q[1];
rz(-1.78396) q[1];
sx q[1];
rz(2.7480385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2834445) q[0];
sx q[0];
rz(-2.1513916) q[0];
sx q[0];
rz(-0.20943187) q[0];
rz(0.15109328) q[2];
sx q[2];
rz(-2.4799621) q[2];
sx q[2];
rz(0.90666319) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4992884) q[1];
sx q[1];
rz(-0.44407108) q[1];
sx q[1];
rz(-1.3127692) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0575074) q[3];
sx q[3];
rz(-1.2197139) q[3];
sx q[3];
rz(0.55817348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5081818) q[2];
sx q[2];
rz(-1.1339302) q[2];
sx q[2];
rz(0.22937648) q[2];
rz(-3.0502099) q[3];
sx q[3];
rz(-1.6978076) q[3];
sx q[3];
rz(2.5858333) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2185739) q[0];
sx q[0];
rz(-0.78892437) q[0];
sx q[0];
rz(-2.5564585) q[0];
rz(1.2773889) q[1];
sx q[1];
rz(-1.2150512) q[1];
sx q[1];
rz(2.9871984) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0421678) q[0];
sx q[0];
rz(-1.4615566) q[0];
sx q[0];
rz(-0.43719284) q[0];
x q[1];
rz(-1.9188668) q[2];
sx q[2];
rz(-1.9744951) q[2];
sx q[2];
rz(0.26223809) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1484376) q[1];
sx q[1];
rz(-1.2497836) q[1];
sx q[1];
rz(1.2953207) q[1];
rz(1.9089237) q[3];
sx q[3];
rz(-2.6498389) q[3];
sx q[3];
rz(-1.7724747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.020891) q[2];
sx q[2];
rz(-0.88940826) q[2];
sx q[2];
rz(-2.6013539) q[2];
rz(0.39515105) q[3];
sx q[3];
rz(-0.65427798) q[3];
sx q[3];
rz(1.7925709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9304792) q[0];
sx q[0];
rz(-1.8607288) q[0];
sx q[0];
rz(1.6092009) q[0];
rz(-2.2156175) q[1];
sx q[1];
rz(-1.7151058) q[1];
sx q[1];
rz(-1.6645974) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84718448) q[0];
sx q[0];
rz(-1.5715181) q[0];
sx q[0];
rz(-0.72615726) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29472827) q[2];
sx q[2];
rz(-1.9325614) q[2];
sx q[2];
rz(-0.28361646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7565347) q[1];
sx q[1];
rz(-1.9137662) q[1];
sx q[1];
rz(0.98481957) q[1];
rz(-0.9524882) q[3];
sx q[3];
rz(-2.455061) q[3];
sx q[3];
rz(1.8812979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8605139) q[2];
sx q[2];
rz(-1.3382567) q[2];
sx q[2];
rz(-2.4930084) q[2];
rz(-0.72566882) q[3];
sx q[3];
rz(-0.23894335) q[3];
sx q[3];
rz(-1.3919977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-2.9072472) q[2];
sx q[2];
rz(-0.6498944) q[2];
sx q[2];
rz(1.9882974) q[2];
rz(0.81448803) q[3];
sx q[3];
rz(-1.9662439) q[3];
sx q[3];
rz(1.6821485) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

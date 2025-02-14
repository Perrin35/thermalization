OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1281066) q[0];
sx q[0];
rz(3.9689316) q[0];
sx q[0];
rz(7.6710424) q[0];
rz(0.85343051) q[1];
sx q[1];
rz(-0.4018468) q[1];
sx q[1];
rz(2.0274577) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3281109) q[0];
sx q[0];
rz(-0.8719647) q[0];
sx q[0];
rz(-0.33150406) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70504712) q[2];
sx q[2];
rz(-1.4234666) q[2];
sx q[2];
rz(1.4718057) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4526032) q[1];
sx q[1];
rz(-1.5509948) q[1];
sx q[1];
rz(0.0048586998) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5680892) q[3];
sx q[3];
rz(-0.57799229) q[3];
sx q[3];
rz(-1.0740394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9806597) q[2];
sx q[2];
rz(-2.7200343) q[2];
sx q[2];
rz(-1.4483615) q[2];
rz(-0.24250044) q[3];
sx q[3];
rz(-0.91553965) q[3];
sx q[3];
rz(-1.2601132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5777957) q[0];
sx q[0];
rz(-1.3160492) q[0];
sx q[0];
rz(-1.0003566) q[0];
rz(0.73161221) q[1];
sx q[1];
rz(-1.7121406) q[1];
sx q[1];
rz(1.9614722) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48413101) q[0];
sx q[0];
rz(-0.13741446) q[0];
sx q[0];
rz(2.4262587) q[0];
x q[1];
rz(1.2928455) q[2];
sx q[2];
rz(-1.5382082) q[2];
sx q[2];
rz(-2.0253092) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0457256) q[1];
sx q[1];
rz(-1.57267) q[1];
sx q[1];
rz(-1.7049403) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3445156) q[3];
sx q[3];
rz(-0.38603544) q[3];
sx q[3];
rz(-2.4222056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0483027) q[2];
sx q[2];
rz(-0.90128171) q[2];
sx q[2];
rz(1.0789336) q[2];
rz(-3.0537187) q[3];
sx q[3];
rz(-1.7885957) q[3];
sx q[3];
rz(-2.8972076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74697772) q[0];
sx q[0];
rz(-1.1629539) q[0];
sx q[0];
rz(-1.8633457) q[0];
rz(-2.6784015) q[1];
sx q[1];
rz(-1.7852424) q[1];
sx q[1];
rz(0.82122222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44959497) q[0];
sx q[0];
rz(-2.2173081) q[0];
sx q[0];
rz(0.23289911) q[0];
rz(-pi) q[1];
rz(-1.5089919) q[2];
sx q[2];
rz(-1.2001749) q[2];
sx q[2];
rz(-0.48423094) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1248904) q[1];
sx q[1];
rz(-1.115106) q[1];
sx q[1];
rz(1.4067006) q[1];
rz(-pi) q[2];
rz(-2.3825112) q[3];
sx q[3];
rz(-2.0584752) q[3];
sx q[3];
rz(1.9239192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.897573) q[2];
sx q[2];
rz(-1.579318) q[2];
sx q[2];
rz(-2.139034) q[2];
rz(-1.2282486) q[3];
sx q[3];
rz(-1.7398261) q[3];
sx q[3];
rz(1.7206515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2380075) q[0];
sx q[0];
rz(-1.0500195) q[0];
sx q[0];
rz(-0.28582698) q[0];
rz(-2.8803275) q[1];
sx q[1];
rz(-1.3328054) q[1];
sx q[1];
rz(3.1108943) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945902) q[0];
sx q[0];
rz(-1.497921) q[0];
sx q[0];
rz(2.8776309) q[0];
rz(0.89392406) q[2];
sx q[2];
rz(-1.6325762) q[2];
sx q[2];
rz(3.1368352) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0847807) q[1];
sx q[1];
rz(-2.0151981) q[1];
sx q[1];
rz(-0.2803352) q[1];
rz(-pi) q[2];
rz(-3.0604393) q[3];
sx q[3];
rz(-1.8834682) q[3];
sx q[3];
rz(-2.816538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3280481) q[2];
sx q[2];
rz(-1.8009461) q[2];
sx q[2];
rz(0.0030585232) q[2];
rz(0.84117738) q[3];
sx q[3];
rz(-2.5523461) q[3];
sx q[3];
rz(-0.35471788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67741126) q[0];
sx q[0];
rz(-0.84742904) q[0];
sx q[0];
rz(1.4671951) q[0];
rz(1.6293619) q[1];
sx q[1];
rz(-2.1758175) q[1];
sx q[1];
rz(-0.95219749) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5291572) q[0];
sx q[0];
rz(-1.5780492) q[0];
sx q[0];
rz(3.0713505) q[0];
rz(-pi) q[1];
rz(1.9604076) q[2];
sx q[2];
rz(-2.2305373) q[2];
sx q[2];
rz(-0.3241186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0972916) q[1];
sx q[1];
rz(-0.11137577) q[1];
sx q[1];
rz(-1.9457413) q[1];
rz(-pi) q[2];
rz(2.2347225) q[3];
sx q[3];
rz(-2.0603486) q[3];
sx q[3];
rz(3.1328894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.78990045) q[2];
sx q[2];
rz(-1.60195) q[2];
sx q[2];
rz(0.44463012) q[2];
rz(2.6897258) q[3];
sx q[3];
rz(-0.63932747) q[3];
sx q[3];
rz(3.0795081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43279466) q[0];
sx q[0];
rz(-2.4169156) q[0];
sx q[0];
rz(-2.2213347) q[0];
rz(-0.96011773) q[1];
sx q[1];
rz(-1.7476387) q[1];
sx q[1];
rz(2.6422909) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2101406) q[0];
sx q[0];
rz(-0.37100756) q[0];
sx q[0];
rz(0.4274344) q[0];
rz(2.9174706) q[2];
sx q[2];
rz(-0.67356985) q[2];
sx q[2];
rz(-2.3526255) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0671406) q[1];
sx q[1];
rz(-2.4765091) q[1];
sx q[1];
rz(2.6159899) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7654789) q[3];
sx q[3];
rz(-0.77738471) q[3];
sx q[3];
rz(0.57143593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.7509191) q[2];
sx q[2];
rz(-1.6030739) q[2];
sx q[2];
rz(-2.1092559) q[2];
rz(0.8647626) q[3];
sx q[3];
rz(-1.7059749) q[3];
sx q[3];
rz(2.5800956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.432935) q[0];
sx q[0];
rz(-1.189804) q[0];
sx q[0];
rz(1.3849965) q[0];
rz(-2.2191091) q[1];
sx q[1];
rz(-1.9360767) q[1];
sx q[1];
rz(0.14499697) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3475048) q[0];
sx q[0];
rz(-1.4055168) q[0];
sx q[0];
rz(1.3896304) q[0];
x q[1];
rz(-0.85741557) q[2];
sx q[2];
rz(-1.5930297) q[2];
sx q[2];
rz(0.93728256) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.26321) q[1];
sx q[1];
rz(-3.0621798) q[1];
sx q[1];
rz(-1.9694049) q[1];
x q[2];
rz(0.77665601) q[3];
sx q[3];
rz(-1.1478416) q[3];
sx q[3];
rz(-2.2222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9292235) q[2];
sx q[2];
rz(-0.23304686) q[2];
sx q[2];
rz(1.2082427) q[2];
rz(-0.82331795) q[3];
sx q[3];
rz(-1.9602785) q[3];
sx q[3];
rz(-1.3594782) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061148297) q[0];
sx q[0];
rz(-2.3539982) q[0];
sx q[0];
rz(1.0795235) q[0];
rz(-0.98006717) q[1];
sx q[1];
rz(-1.020224) q[1];
sx q[1];
rz(-0.10990873) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82775527) q[0];
sx q[0];
rz(-1.8640564) q[0];
sx q[0];
rz(0.37773962) q[0];
x q[1];
rz(-2.8632322) q[2];
sx q[2];
rz(-1.5621857) q[2];
sx q[2];
rz(-1.4573163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3802783) q[1];
sx q[1];
rz(-1.8683234) q[1];
sx q[1];
rz(2.6975836) q[1];
x q[2];
rz(2.1004543) q[3];
sx q[3];
rz(-1.7653437) q[3];
sx q[3];
rz(1.1916849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2073233) q[2];
sx q[2];
rz(-2.0685652) q[2];
sx q[2];
rz(2.4369924) q[2];
rz(-2.8568824) q[3];
sx q[3];
rz(-2.4094818) q[3];
sx q[3];
rz(2.708191) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1339742) q[0];
sx q[0];
rz(-1.3173137) q[0];
sx q[0];
rz(-2.2229069) q[0];
rz(-0.48565117) q[1];
sx q[1];
rz(-2.698027) q[1];
sx q[1];
rz(-1.1007016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9022517) q[0];
sx q[0];
rz(-1.709769) q[0];
sx q[0];
rz(-0.13987565) q[0];
rz(-pi) q[1];
rz(1.2687911) q[2];
sx q[2];
rz(-1.3493611) q[2];
sx q[2];
rz(2.1083652) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6599258) q[1];
sx q[1];
rz(-1.9710014) q[1];
sx q[1];
rz(-0.11379644) q[1];
rz(-pi) q[2];
rz(1.507255) q[3];
sx q[3];
rz(-2.2054234) q[3];
sx q[3];
rz(1.612048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47621581) q[2];
sx q[2];
rz(-2.2874338) q[2];
sx q[2];
rz(0.7594792) q[2];
rz(-2.1935943) q[3];
sx q[3];
rz(-1.6839323) q[3];
sx q[3];
rz(1.6215526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(3.0976902) q[0];
sx q[0];
rz(-0.49711415) q[0];
sx q[0];
rz(-2.8105766) q[0];
rz(2.379592) q[1];
sx q[1];
rz(-0.29251978) q[1];
sx q[1];
rz(2.5133572) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8593438) q[0];
sx q[0];
rz(-2.4971136) q[0];
sx q[0];
rz(-0.47596875) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8339001) q[2];
sx q[2];
rz(-1.7950222) q[2];
sx q[2];
rz(-2.4196408) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8208295) q[1];
sx q[1];
rz(-0.71457982) q[1];
sx q[1];
rz(0.54820593) q[1];
rz(-pi) q[2];
rz(0.41251934) q[3];
sx q[3];
rz(-1.8897353) q[3];
sx q[3];
rz(2.7273415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1198173) q[2];
sx q[2];
rz(-0.47406083) q[2];
sx q[2];
rz(-0.21381703) q[2];
rz(2.1553433) q[3];
sx q[3];
rz(-1.2912368) q[3];
sx q[3];
rz(0.10828644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66961359) q[0];
sx q[0];
rz(-1.6983953) q[0];
sx q[0];
rz(-1.4629913) q[0];
rz(-1.1769453) q[1];
sx q[1];
rz(-1.4851478) q[1];
sx q[1];
rz(-3.1335395) q[1];
rz(1.2659579) q[2];
sx q[2];
rz(-1.4043456) q[2];
sx q[2];
rz(3.059491) q[2];
rz(0.65153788) q[3];
sx q[3];
rz(-1.6350411) q[3];
sx q[3];
rz(-1.7512334) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

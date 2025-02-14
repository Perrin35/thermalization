OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2562113) q[0];
sx q[0];
rz(5.5698759) q[0];
sx q[0];
rz(8.2110693) q[0];
rz(-1.3656536) q[1];
sx q[1];
rz(-2.8627099) q[1];
sx q[1];
rz(-0.20067781) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532625) q[0];
sx q[0];
rz(-1.6536667) q[0];
sx q[0];
rz(-1.1088662) q[0];
x q[1];
rz(-0.63849475) q[2];
sx q[2];
rz(-0.26836499) q[2];
sx q[2];
rz(0.61682781) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7534448) q[1];
sx q[1];
rz(-0.26784836) q[1];
sx q[1];
rz(2.4004188) q[1];
x q[2];
rz(2.6418125) q[3];
sx q[3];
rz(-0.3480787) q[3];
sx q[3];
rz(-1.6655371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.063244907) q[2];
sx q[2];
rz(-1.1996317) q[2];
sx q[2];
rz(2.4285748) q[2];
rz(1.2837422) q[3];
sx q[3];
rz(-2.4132437) q[3];
sx q[3];
rz(2.242105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9707608) q[0];
sx q[0];
rz(-1.66667) q[0];
sx q[0];
rz(-1.4738039) q[0];
rz(-2.5693192) q[1];
sx q[1];
rz(-2.0794561) q[1];
sx q[1];
rz(-1.3776113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.361918) q[0];
sx q[0];
rz(-2.5847844) q[0];
sx q[0];
rz(-0.98018719) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6102561) q[2];
sx q[2];
rz(-1.1684993) q[2];
sx q[2];
rz(0.50237331) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26220825) q[1];
sx q[1];
rz(-1.9010494) q[1];
sx q[1];
rz(0.31008215) q[1];
rz(-pi) q[2];
rz(0.20444025) q[3];
sx q[3];
rz(-1.8722765) q[3];
sx q[3];
rz(0.35273509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1231692) q[2];
sx q[2];
rz(-1.6320684) q[2];
sx q[2];
rz(1.8040166) q[2];
rz(0.34559524) q[3];
sx q[3];
rz(-0.59185043) q[3];
sx q[3];
rz(-0.93446294) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0551986) q[0];
sx q[0];
rz(-0.18562695) q[0];
sx q[0];
rz(0.58919543) q[0];
rz(-0.20801726) q[1];
sx q[1];
rz(-2.5841525) q[1];
sx q[1];
rz(2.9357108) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.170144) q[0];
sx q[0];
rz(-1.2177694) q[0];
sx q[0];
rz(2.5021509) q[0];
rz(-0.85742204) q[2];
sx q[2];
rz(-1.1829585) q[2];
sx q[2];
rz(2.037354) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.65218788) q[1];
sx q[1];
rz(-0.88880782) q[1];
sx q[1];
rz(-1.1754065) q[1];
x q[2];
rz(0.08783665) q[3];
sx q[3];
rz(-1.1655775) q[3];
sx q[3];
rz(-2.3072568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0516633) q[2];
sx q[2];
rz(-1.0445003) q[2];
sx q[2];
rz(2.5970686) q[2];
rz(-0.45474592) q[3];
sx q[3];
rz(-2.5179722) q[3];
sx q[3];
rz(3.1375942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74314463) q[0];
sx q[0];
rz(-0.54773098) q[0];
sx q[0];
rz(2.4416583) q[0];
rz(-1.3088538) q[1];
sx q[1];
rz(-1.0700285) q[1];
sx q[1];
rz(1.4189789) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0467906) q[0];
sx q[0];
rz(-1.3700458) q[0];
sx q[0];
rz(-1.0958452) q[0];
rz(-pi) q[1];
rz(-0.019446418) q[2];
sx q[2];
rz(-2.0692973) q[2];
sx q[2];
rz(0.74048238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5955453) q[1];
sx q[1];
rz(-0.7460119) q[1];
sx q[1];
rz(2.7849848) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1484234) q[3];
sx q[3];
rz(-2.7971533) q[3];
sx q[3];
rz(0.41448739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2122638) q[2];
sx q[2];
rz(-0.45390359) q[2];
sx q[2];
rz(1.7960499) q[2];
rz(-2.4228607) q[3];
sx q[3];
rz(-2.4001382) q[3];
sx q[3];
rz(-0.68812266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5162002) q[0];
sx q[0];
rz(-1.9560408) q[0];
sx q[0];
rz(-2.5984883) q[0];
rz(-2.886046) q[1];
sx q[1];
rz(-2.6592022) q[1];
sx q[1];
rz(1.7593174) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10777625) q[0];
sx q[0];
rz(-1.9703426) q[0];
sx q[0];
rz(-2.5679213) q[0];
x q[1];
rz(-0.15145258) q[2];
sx q[2];
rz(-0.77255682) q[2];
sx q[2];
rz(0.59237769) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71798793) q[1];
sx q[1];
rz(-1.8306872) q[1];
sx q[1];
rz(2.1433448) q[1];
rz(-pi) q[2];
rz(2.7543254) q[3];
sx q[3];
rz(-2.7125689) q[3];
sx q[3];
rz(2.9128592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9690659) q[2];
sx q[2];
rz(-1.669599) q[2];
sx q[2];
rz(0.36953163) q[2];
rz(-1.8116123) q[3];
sx q[3];
rz(-2.3510635) q[3];
sx q[3];
rz(1.3548405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1177649) q[0];
sx q[0];
rz(-0.33482877) q[0];
sx q[0];
rz(-0.086300015) q[0];
rz(1.49508) q[1];
sx q[1];
rz(-1.9638655) q[1];
sx q[1];
rz(-2.7118447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54329007) q[0];
sx q[0];
rz(-1.5006646) q[0];
sx q[0];
rz(-1.381641) q[0];
rz(-pi) q[1];
rz(0.65799539) q[2];
sx q[2];
rz(-1.3086196) q[2];
sx q[2];
rz(0.10181759) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8064902) q[1];
sx q[1];
rz(-1.1336328) q[1];
sx q[1];
rz(-2.6173068) q[1];
rz(-0.67062258) q[3];
sx q[3];
rz(-1.3486996) q[3];
sx q[3];
rz(0.64034739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0325844) q[2];
sx q[2];
rz(-1.0823559) q[2];
sx q[2];
rz(1.2004131) q[2];
rz(-0.16610185) q[3];
sx q[3];
rz(-2.3717272) q[3];
sx q[3];
rz(-2.2782245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2849041) q[0];
sx q[0];
rz(-2.3785474) q[0];
sx q[0];
rz(-0.82409182) q[0];
rz(1.2245945) q[1];
sx q[1];
rz(-1.2242182) q[1];
sx q[1];
rz(-2.1276988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4201502) q[0];
sx q[0];
rz(-2.3334896) q[0];
sx q[0];
rz(1.1366913) q[0];
rz(-pi) q[1];
rz(-2.7043506) q[2];
sx q[2];
rz(-1.8595942) q[2];
sx q[2];
rz(-1.0693969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6142004) q[1];
sx q[1];
rz(-1.0353312) q[1];
sx q[1];
rz(2.0245069) q[1];
x q[2];
rz(2.7498662) q[3];
sx q[3];
rz(-2.0731032) q[3];
sx q[3];
rz(1.6592178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8635204) q[2];
sx q[2];
rz(-0.72276989) q[2];
sx q[2];
rz(-1.8966759) q[2];
rz(-0.23166367) q[3];
sx q[3];
rz(-1.04117) q[3];
sx q[3];
rz(2.9932573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046655) q[0];
sx q[0];
rz(-1.2060839) q[0];
sx q[0];
rz(-0.99223247) q[0];
rz(-0.16920432) q[1];
sx q[1];
rz(-2.2997586) q[1];
sx q[1];
rz(-1.7281035) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.468576) q[0];
sx q[0];
rz(-1.6454738) q[0];
sx q[0];
rz(-0.59665307) q[0];
rz(-pi) q[1];
rz(1.1216759) q[2];
sx q[2];
rz(-2.5418021) q[2];
sx q[2];
rz(3.1315294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1991829) q[1];
sx q[1];
rz(-1.8687667) q[1];
sx q[1];
rz(-1.6040785) q[1];
x q[2];
rz(-2.0510717) q[3];
sx q[3];
rz(-1.3930219) q[3];
sx q[3];
rz(-0.67713523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8168872) q[2];
sx q[2];
rz(-1.6605261) q[2];
sx q[2];
rz(1.4229577) q[2];
rz(-0.11296806) q[3];
sx q[3];
rz(-0.26765099) q[3];
sx q[3];
rz(2.512219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9519575) q[0];
sx q[0];
rz(-1.4653787) q[0];
sx q[0];
rz(-2.5919609) q[0];
rz(2.5426087) q[1];
sx q[1];
rz(-0.87268972) q[1];
sx q[1];
rz(-1.6302861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0881168) q[0];
sx q[0];
rz(-1.5673182) q[0];
sx q[0];
rz(2.9649078) q[0];
rz(-pi) q[1];
rz(-2.6312066) q[2];
sx q[2];
rz(-0.56737075) q[2];
sx q[2];
rz(-2.9387143) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.05311) q[1];
sx q[1];
rz(-2.4743926) q[1];
sx q[1];
rz(-0.88796292) q[1];
x q[2];
rz(-0.49137791) q[3];
sx q[3];
rz(-1.8432321) q[3];
sx q[3];
rz(0.85173729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0281684) q[2];
sx q[2];
rz(-1.3415965) q[2];
sx q[2];
rz(2.8070731) q[2];
rz(-0.64535514) q[3];
sx q[3];
rz(-2.1367475) q[3];
sx q[3];
rz(2.9042802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065011218) q[0];
sx q[0];
rz(-2.8028221) q[0];
sx q[0];
rz(2.5184799) q[0];
rz(-2.8399641) q[1];
sx q[1];
rz(-2.5239065) q[1];
sx q[1];
rz(-1.041144) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89462507) q[0];
sx q[0];
rz(-2.0026221) q[0];
sx q[0];
rz(0.65460848) q[0];
rz(-0.1847965) q[2];
sx q[2];
rz(-0.46773887) q[2];
sx q[2];
rz(-0.50450215) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.54339441) q[1];
sx q[1];
rz(-2.3375727) q[1];
sx q[1];
rz(-1.5280185) q[1];
rz(-pi) q[2];
rz(-2.9277293) q[3];
sx q[3];
rz(-1.2215316) q[3];
sx q[3];
rz(1.6598778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4447896) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(-0.94318548) q[2];
rz(1.7275564) q[3];
sx q[3];
rz(-0.90641886) q[3];
sx q[3];
rz(2.3448155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4010451) q[0];
sx q[0];
rz(-2.7475806) q[0];
sx q[0];
rz(-0.077234118) q[0];
rz(2.7229591) q[1];
sx q[1];
rz(-1.5593465) q[1];
sx q[1];
rz(-2.9895463) q[1];
rz(0.05337333) q[2];
sx q[2];
rz(-0.57562258) q[2];
sx q[2];
rz(-2.876566) q[2];
rz(-0.54316212) q[3];
sx q[3];
rz(-0.49351963) q[3];
sx q[3];
rz(-1.7388572) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

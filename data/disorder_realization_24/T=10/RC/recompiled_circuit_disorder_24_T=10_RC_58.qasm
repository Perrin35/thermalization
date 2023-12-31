OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(-1.0868602) q[0];
sx q[0];
rz(-1.342919) q[0];
rz(7.8328447) q[1];
sx q[1];
rz(3.0631493) q[1];
sx q[1];
rz(6.9258239) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6519878) q[0];
sx q[0];
rz(-1.542924) q[0];
sx q[0];
rz(-0.53919381) q[0];
rz(-2.6334555) q[2];
sx q[2];
rz(-1.7388441) q[2];
sx q[2];
rz(-0.44067581) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.77036422) q[1];
sx q[1];
rz(-0.1703573) q[1];
sx q[1];
rz(0.78070663) q[1];
rz(-pi) q[2];
rz(1.5343127) q[3];
sx q[3];
rz(-1.8047223) q[3];
sx q[3];
rz(0.28947383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(-1.2791963) q[2];
rz(0.71875087) q[3];
sx q[3];
rz(-1.5703392) q[3];
sx q[3];
rz(2.8180502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779125) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(-0.15788831) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(0.78871361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65608998) q[0];
sx q[0];
rz(-2.9733109) q[0];
sx q[0];
rz(2.3165354) q[0];
rz(3.0208203) q[2];
sx q[2];
rz(-2.883652) q[2];
sx q[2];
rz(-2.068012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5309108) q[1];
sx q[1];
rz(-1.0265988) q[1];
sx q[1];
rz(1.9287841) q[1];
rz(-pi) q[2];
rz(-2.4103568) q[3];
sx q[3];
rz(-0.9160708) q[3];
sx q[3];
rz(0.54192858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7754037) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(-0.186084) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(-0.11793605) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9300951) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(0.63013664) q[0];
rz(3.0139626) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(0.72174597) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16746178) q[0];
sx q[0];
rz(-0.92739096) q[0];
sx q[0];
rz(1.2533623) q[0];
rz(-pi) q[1];
rz(2.6038405) q[2];
sx q[2];
rz(-0.76965145) q[2];
sx q[2];
rz(-0.69307454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7194781) q[1];
sx q[1];
rz(-2.0146703) q[1];
sx q[1];
rz(1.5281954) q[1];
x q[2];
rz(2.0201163) q[3];
sx q[3];
rz(-2.46393) q[3];
sx q[3];
rz(0.47798702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7599941) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-2.2325113) q[2];
rz(-0.51820731) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58406126) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(2.361239) q[1];
sx q[1];
rz(-0.50023729) q[1];
sx q[1];
rz(-0.76400486) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4064179) q[0];
sx q[0];
rz(-1.5582726) q[0];
sx q[0];
rz(-2.2494715) q[0];
x q[1];
rz(-2.4970384) q[2];
sx q[2];
rz(-1.6719712) q[2];
sx q[2];
rz(-0.25601706) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3180247) q[1];
sx q[1];
rz(-1.9293702) q[1];
sx q[1];
rz(-2.0477247) q[1];
rz(-0.56960168) q[3];
sx q[3];
rz(-0.74865018) q[3];
sx q[3];
rz(-0.99018712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8445231) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(2.6085473) q[2];
rz(-0.26238966) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2247291) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(1.2438783) q[0];
rz(-0.22661701) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(2.7382543) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82542244) q[0];
sx q[0];
rz(-0.24992019) q[0];
sx q[0];
rz(0.73096801) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7062347) q[2];
sx q[2];
rz(-2.407496) q[2];
sx q[2];
rz(-1.8045319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7008608) q[1];
sx q[1];
rz(-1.7248389) q[1];
sx q[1];
rz(0.53149077) q[1];
rz(-pi) q[2];
rz(-1.0395398) q[3];
sx q[3];
rz(-2.1752393) q[3];
sx q[3];
rz(-1.2951375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.32101813) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(2.3590951) q[2];
rz(2.0292422) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.430442) q[0];
sx q[0];
rz(-2.7395881) q[0];
sx q[0];
rz(2.6859786) q[0];
rz(-3.0420711) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(-3.1351556) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59505263) q[0];
sx q[0];
rz(-1.8227302) q[0];
sx q[0];
rz(-0.76058723) q[0];
x q[1];
rz(1.3077277) q[2];
sx q[2];
rz(-0.38882133) q[2];
sx q[2];
rz(0.013465492) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.41758075) q[1];
sx q[1];
rz(-1.9829826) q[1];
sx q[1];
rz(-2.8860303) q[1];
x q[2];
rz(-0.75565831) q[3];
sx q[3];
rz(-2.14058) q[3];
sx q[3];
rz(0.33982402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(0.84645611) q[2];
rz(-0.99772292) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098175123) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(0.78272351) q[1];
sx q[1];
rz(-2.0109773) q[1];
sx q[1];
rz(1.9810716) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6877277) q[0];
sx q[0];
rz(-1.9181607) q[0];
sx q[0];
rz(1.0930644) q[0];
rz(0.47126982) q[2];
sx q[2];
rz(-1.0980714) q[2];
sx q[2];
rz(-1.3408692) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6351663) q[1];
sx q[1];
rz(-0.70235683) q[1];
sx q[1];
rz(0.96794767) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7619744) q[3];
sx q[3];
rz(-1.2109204) q[3];
sx q[3];
rz(0.66756638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(-0.78545061) q[2];
rz(-2.3857332) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(-2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3128368) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(1.998741) q[0];
rz(1.3061334) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(-0.41608861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0932255) q[0];
sx q[0];
rz(-1.0972293) q[0];
sx q[0];
rz(2.918539) q[0];
rz(-pi) q[1];
rz(2.9990254) q[2];
sx q[2];
rz(-1.1066184) q[2];
sx q[2];
rz(-0.67205059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1382653) q[1];
sx q[1];
rz(-1.4107553) q[1];
sx q[1];
rz(0.40469594) q[1];
rz(-pi) q[2];
rz(-3.0307426) q[3];
sx q[3];
rz(-1.1006315) q[3];
sx q[3];
rz(1.6925616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0351506) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(-0.32315928) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(-2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(1.9966104) q[0];
rz(-1.9454983) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(-2.6224565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0275092) q[0];
sx q[0];
rz(-0.4013831) q[0];
sx q[0];
rz(-1.2252349) q[0];
x q[1];
rz(0.81705117) q[2];
sx q[2];
rz(-2.8231986) q[2];
sx q[2];
rz(1.5621834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0172826) q[1];
sx q[1];
rz(-0.78622422) q[1];
sx q[1];
rz(-0.12676523) q[1];
rz(-pi) q[2];
rz(-1.395123) q[3];
sx q[3];
rz(-1.3364949) q[3];
sx q[3];
rz(-2.7845886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8273932) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(-3.1006151) q[2];
rz(2.273902) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7460019) q[0];
sx q[0];
rz(-2.262291) q[0];
sx q[0];
rz(1.6145153) q[0];
rz(1.4279667) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(1.6428927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258543) q[0];
sx q[0];
rz(-0.07809528) q[0];
sx q[0];
rz(1.8321091) q[0];
x q[1];
rz(-1.2730359) q[2];
sx q[2];
rz(-0.16212633) q[2];
sx q[2];
rz(-0.33725421) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5012706) q[1];
sx q[1];
rz(-2.3295799) q[1];
sx q[1];
rz(-1.0971607) q[1];
rz(-pi) q[2];
rz(2.5999971) q[3];
sx q[3];
rz(-0.23532) q[3];
sx q[3];
rz(-3.1129587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3830118) q[2];
sx q[2];
rz(-2.0643533) q[2];
sx q[2];
rz(-0.14653462) q[2];
rz(-2.3274029) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(2.7248785) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9183337) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(1.2659484) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(-2.0064034) q[2];
sx q[2];
rz(-0.26293593) q[2];
sx q[2];
rz(1.5756366) q[2];
rz(2.2371348) q[3];
sx q[3];
rz(-0.96257985) q[3];
sx q[3];
rz(-1.0812159) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

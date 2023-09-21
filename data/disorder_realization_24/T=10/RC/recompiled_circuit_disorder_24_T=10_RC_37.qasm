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
rz(2.0547325) q[0];
sx q[0];
rz(7.6261043) q[0];
rz(7.8328447) q[1];
sx q[1];
rz(3.0631493) q[1];
sx q[1];
rz(6.9258239) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064518236) q[0];
sx q[0];
rz(-1.0318349) q[0];
sx q[0];
rz(-1.6032739) q[0];
rz(-2.8060993) q[2];
sx q[2];
rz(-0.53288424) q[2];
sx q[2];
rz(2.3032308) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77036422) q[1];
sx q[1];
rz(-0.1703573) q[1];
sx q[1];
rz(-0.78070663) q[1];
rz(-pi) q[2];
x q[2];
rz(1.60728) q[3];
sx q[3];
rz(-1.8047223) q[3];
sx q[3];
rz(-0.28947383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(-1.8623964) q[2];
rz(-0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(2.8180502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779125) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(-0.15788831) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(2.352879) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6533587) q[0];
sx q[0];
rz(-1.6846859) q[0];
sx q[0];
rz(1.4466404) q[0];
x q[1];
rz(-2.8854495) q[2];
sx q[2];
rz(-1.5400585) q[2];
sx q[2];
rz(2.5275633) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1566369) q[1];
sx q[1];
rz(-2.5002694) q[1];
sx q[1];
rz(-2.6167469) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85487811) q[3];
sx q[3];
rz(-0.93920556) q[3];
sx q[3];
rz(-2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7754037) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(-2.9555087) q[2];
rz(-0.6535334) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2114975) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(-2.511456) q[0];
rz(0.12763003) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(0.72174597) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9741309) q[0];
sx q[0];
rz(-0.92739096) q[0];
sx q[0];
rz(1.2533623) q[0];
x q[1];
rz(-1.1100936) q[2];
sx q[2];
rz(-0.93020541) q[2];
sx q[2];
rz(-5.0355807e-05) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.1303781) q[1];
sx q[1];
rz(-1.6092669) q[1];
sx q[1];
rz(-0.4442261) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8053022) q[3];
sx q[3];
rz(-2.170917) q[3];
sx q[3];
rz(1.0322514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3815986) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-2.2325113) q[2];
rz(2.6233853) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(-0.28809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5575314) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(-2.361239) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(2.3775878) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873136) q[0];
sx q[0];
rz(-2.2494082) q[0];
sx q[0];
rz(-0.016088386) q[0];
rz(-pi) q[1];
rz(-1.4444703) q[2];
sx q[2];
rz(-2.2115123) q[2];
sx q[2];
rz(1.7510406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82356794) q[1];
sx q[1];
rz(-1.2122224) q[1];
sx q[1];
rz(-1.093868) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4776811) q[3];
sx q[3];
rz(-1.1949364) q[3];
sx q[3];
rz(1.0192878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(-0.53304535) q[2];
rz(-0.26238966) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9168636) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(1.8977144) q[0];
rz(2.9149756) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(0.40333834) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3161702) q[0];
sx q[0];
rz(-2.8916725) q[0];
sx q[0];
rz(2.4106246) q[0];
rz(-pi) q[1];
rz(0.43535797) q[2];
sx q[2];
rz(-2.407496) q[2];
sx q[2];
rz(1.3370607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9215645) q[1];
sx q[1];
rz(-1.0462531) q[1];
sx q[1];
rz(1.7490053) q[1];
x q[2];
rz(-0.63286085) q[3];
sx q[3];
rz(-0.78213464) q[3];
sx q[3];
rz(0.49304214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32101813) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(-2.3590951) q[2];
rz(-2.0292422) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(-2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(2.6859786) q[0];
rz(-3.0420711) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(3.1351556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71953668) q[0];
sx q[0];
rz(-0.7932084) q[0];
sx q[0];
rz(2.7842245) q[0];
rz(-3.0354584) q[2];
sx q[2];
rz(-1.1960293) q[2];
sx q[2];
rz(-2.8449164) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7240119) q[1];
sx q[1];
rz(-1.9829826) q[1];
sx q[1];
rz(0.25556232) q[1];
x q[2];
rz(2.3859343) q[3];
sx q[3];
rz(-2.14058) q[3];
sx q[3];
rz(0.33982402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36879888) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(-2.2951365) q[2];
rz(2.1438697) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0434175) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(-0.78272351) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.9810716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6877277) q[0];
sx q[0];
rz(-1.9181607) q[0];
sx q[0];
rz(-2.0485282) q[0];
rz(-pi) q[1];
rz(2.6703228) q[2];
sx q[2];
rz(-1.0980714) q[2];
sx q[2];
rz(-1.8007235) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5481422) q[1];
sx q[1];
rz(-1.1957809) q[1];
sx q[1];
rz(2.1795991) q[1];
rz(-pi) q[2];
rz(1.9556932) q[3];
sx q[3];
rz(-1.2166096) q[3];
sx q[3];
rz(2.3779496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.065585) q[2];
sx q[2];
rz(-0.92131725) q[2];
sx q[2];
rz(2.356142) q[2];
rz(-2.3857332) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(0.41539645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3128368) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(-1.998741) q[0];
rz(-1.3061334) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(-0.41608861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4131665) q[0];
sx q[0];
rz(-2.6217555) q[0];
sx q[0];
rz(1.9782938) q[0];
rz(1.2942737) q[2];
sx q[2];
rz(-0.48404901) q[2];
sx q[2];
rz(-2.7798142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63562288) q[1];
sx q[1];
rz(-1.1715679) q[1];
sx q[1];
rz(1.7446306) q[1];
x q[2];
rz(2.0434521) q[3];
sx q[3];
rz(-1.4720159) q[3];
sx q[3];
rz(-0.071382513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1064421) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(-2.9390826) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(0.57730738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(-1.1449822) q[0];
rz(1.1960944) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(0.51913613) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0275092) q[0];
sx q[0];
rz(-2.7402096) q[0];
sx q[0];
rz(-1.9163577) q[0];
rz(1.8066508) q[2];
sx q[2];
rz(-1.354885) q[2];
sx q[2];
rz(2.4049135) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.054159315) q[1];
sx q[1];
rz(-2.3490153) q[1];
sx q[1];
rz(-1.4448318) q[1];
rz(2.9037644) q[3];
sx q[3];
rz(-1.3999709) q[3];
sx q[3];
rz(1.886614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8273932) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(3.1006151) q[2];
rz(-0.86769062) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(-0.51122558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39559078) q[0];
sx q[0];
rz(-2.262291) q[0];
sx q[0];
rz(1.6145153) q[0];
rz(1.7136259) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(1.4987) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.98793) q[0];
sx q[0];
rz(-1.4953574) q[0];
sx q[0];
rz(0.020214202) q[0];
rz(-1.7259049) q[2];
sx q[2];
rz(-1.6181706) q[2];
sx q[2];
rz(0.93946379) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8720819) q[1];
sx q[1];
rz(-1.9081569) q[1];
sx q[1];
rz(-0.81706394) q[1];
x q[2];
rz(-2.9389631) q[3];
sx q[3];
rz(-1.6912795) q[3];
sx q[3];
rz(-1.0128302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3830118) q[2];
sx q[2];
rz(-2.0643533) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(0.81418973) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(-2.7248785) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9183337) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(-1.2659484) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(0.11309359) q[2];
sx q[2];
rz(-1.332915) q[2];
sx q[2];
rz(-1.1168196) q[2];
rz(-0.72487763) q[3];
sx q[3];
rz(-1.0387883) q[3];
sx q[3];
rz(-3.0742857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

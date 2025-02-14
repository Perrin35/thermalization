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
rz(-0.80225575) q[0];
sx q[0];
rz(-1.7576317) q[0];
sx q[0];
rz(-1.8729762) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(-3.0825244) q[1];
sx q[1];
rz(-1.7360092) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48412886) q[0];
sx q[0];
rz(-2.0031558) q[0];
sx q[0];
rz(-2.0922489) q[0];
x q[1];
rz(0.28702847) q[2];
sx q[2];
rz(-0.69106709) q[2];
sx q[2];
rz(-0.075714684) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.51912145) q[1];
sx q[1];
rz(-2.2833385) q[1];
sx q[1];
rz(2.126241) q[1];
rz(-0.22698347) q[3];
sx q[3];
rz(-0.90714083) q[3];
sx q[3];
rz(0.047078156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9872226) q[2];
sx q[2];
rz(-1.1673678) q[2];
sx q[2];
rz(-0.78987375) q[2];
rz(3.1176944) q[3];
sx q[3];
rz(-1.3393211) q[3];
sx q[3];
rz(-0.53643119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.8928878) q[0];
sx q[0];
rz(-1.6510115) q[0];
sx q[0];
rz(1.254068) q[0];
rz(1.4887384) q[1];
sx q[1];
rz(-2.2658927) q[1];
sx q[1];
rz(2.658433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2452908) q[0];
sx q[0];
rz(-1.622513) q[0];
sx q[0];
rz(1.381676) q[0];
x q[1];
rz(0.56098361) q[2];
sx q[2];
rz(-0.81811726) q[2];
sx q[2];
rz(2.4647692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3846783) q[1];
sx q[1];
rz(-1.9325629) q[1];
sx q[1];
rz(-1.2902618) q[1];
x q[2];
rz(-1.7464569) q[3];
sx q[3];
rz(-1.2845728) q[3];
sx q[3];
rz(1.3940879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.88863215) q[2];
sx q[2];
rz(-1.4560207) q[2];
sx q[2];
rz(1.6271094) q[2];
rz(-3.0857981) q[3];
sx q[3];
rz(-1.5443708) q[3];
sx q[3];
rz(-1.6519215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2273939) q[0];
sx q[0];
rz(-2.6970503) q[0];
sx q[0];
rz(3.0134873) q[0];
rz(-2.5654492) q[1];
sx q[1];
rz(-0.73155254) q[1];
sx q[1];
rz(2.3760956) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.573581) q[0];
sx q[0];
rz(-2.4666602) q[0];
sx q[0];
rz(0.57008596) q[0];
rz(1.9472935) q[2];
sx q[2];
rz(-0.91937477) q[2];
sx q[2];
rz(2.5848856) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0772503) q[1];
sx q[1];
rz(-2.0301452) q[1];
sx q[1];
rz(-2.5959754) q[1];
rz(-pi) q[2];
rz(0.14344826) q[3];
sx q[3];
rz(-1.5960403) q[3];
sx q[3];
rz(0.48893602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8219882) q[2];
sx q[2];
rz(-2.2687948) q[2];
sx q[2];
rz(2.2785211) q[2];
rz(2.7989164) q[3];
sx q[3];
rz(-0.42049146) q[3];
sx q[3];
rz(-1.6531403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89219013) q[0];
sx q[0];
rz(-2.1718195) q[0];
sx q[0];
rz(2.6595907) q[0];
rz(0.30872289) q[1];
sx q[1];
rz(-2.8161507) q[1];
sx q[1];
rz(0.6122922) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21633241) q[0];
sx q[0];
rz(-1.060492) q[0];
sx q[0];
rz(0.55522529) q[0];
rz(-pi) q[1];
rz(-3.127549) q[2];
sx q[2];
rz(-1.6545452) q[2];
sx q[2];
rz(1.75911) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0902113) q[1];
sx q[1];
rz(-1.6156229) q[1];
sx q[1];
rz(-0.96446891) q[1];
rz(0.14928603) q[3];
sx q[3];
rz(-1.5642435) q[3];
sx q[3];
rz(0.1803151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65773949) q[2];
sx q[2];
rz(-2.190399) q[2];
sx q[2];
rz(1.8335906) q[2];
rz(-0.43404239) q[3];
sx q[3];
rz(-1.3515892) q[3];
sx q[3];
rz(1.2691809) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16172116) q[0];
sx q[0];
rz(-1.5890108) q[0];
sx q[0];
rz(-0.27444926) q[0];
rz(1.6203923) q[1];
sx q[1];
rz(-1.043964) q[1];
sx q[1];
rz(-1.8843947) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3109821) q[0];
sx q[0];
rz(-0.3468967) q[0];
sx q[0];
rz(-0.78126379) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36894746) q[2];
sx q[2];
rz(-1.0668179) q[2];
sx q[2];
rz(1.1541909) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3538227) q[1];
sx q[1];
rz(-2.1102043) q[1];
sx q[1];
rz(-0.25168519) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2755695) q[3];
sx q[3];
rz(-1.3228205) q[3];
sx q[3];
rz(-3.0803404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3955128) q[2];
sx q[2];
rz(-1.6359676) q[2];
sx q[2];
rz(2.5858322) q[2];
rz(-0.79103509) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(-1.5033495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24285862) q[0];
sx q[0];
rz(-1.4396311) q[0];
sx q[0];
rz(0.71204251) q[0];
rz(2.5391319) q[1];
sx q[1];
rz(-2.7280877) q[1];
sx q[1];
rz(3.0625524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3369737) q[0];
sx q[0];
rz(-0.078402407) q[0];
sx q[0];
rz(0.96167643) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84509648) q[2];
sx q[2];
rz(-0.11339408) q[2];
sx q[2];
rz(1.5528983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1823719) q[1];
sx q[1];
rz(-2.3010265) q[1];
sx q[1];
rz(1.8664136) q[1];
rz(-0.2980663) q[3];
sx q[3];
rz(-1.4722927) q[3];
sx q[3];
rz(2.5677266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.041212335) q[2];
sx q[2];
rz(-0.93116394) q[2];
sx q[2];
rz(0.97990123) q[2];
rz(1.529871) q[3];
sx q[3];
rz(-1.2414705) q[3];
sx q[3];
rz(-2.8821168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3945382) q[0];
sx q[0];
rz(-1.8600445) q[0];
sx q[0];
rz(-0.10093149) q[0];
rz(2.586567) q[1];
sx q[1];
rz(-0.9340159) q[1];
sx q[1];
rz(1.8468044) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053558401) q[0];
sx q[0];
rz(-2.4248666) q[0];
sx q[0];
rz(-1.2117366) q[0];
rz(-0.00097043911) q[2];
sx q[2];
rz(-1.8651267) q[2];
sx q[2];
rz(0.4145588) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.50404149) q[1];
sx q[1];
rz(-0.45256549) q[1];
sx q[1];
rz(2.8082147) q[1];
x q[2];
rz(-1.0247308) q[3];
sx q[3];
rz(-2.6203558) q[3];
sx q[3];
rz(-1.8841528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9635705) q[2];
sx q[2];
rz(-2.2691085) q[2];
sx q[2];
rz(-2.8998609) q[2];
rz(-2.669615) q[3];
sx q[3];
rz(-2.9265672) q[3];
sx q[3];
rz(1.5579461) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80411512) q[0];
sx q[0];
rz(-0.98144704) q[0];
sx q[0];
rz(2.5313983) q[0];
rz(-0.21774165) q[1];
sx q[1];
rz(-2.1861031) q[1];
sx q[1];
rz(2.311923) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10388869) q[0];
sx q[0];
rz(-1.0173326) q[0];
sx q[0];
rz(-2.3174556) q[0];
rz(-pi) q[1];
rz(-1.2118633) q[2];
sx q[2];
rz(-1.8808489) q[2];
sx q[2];
rz(2.884825) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9431994) q[1];
sx q[1];
rz(-1.7560609) q[1];
sx q[1];
rz(2.7279502) q[1];
x q[2];
rz(3.025981) q[3];
sx q[3];
rz(-0.86263166) q[3];
sx q[3];
rz(0.96683217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3024451) q[2];
sx q[2];
rz(-1.382901) q[2];
sx q[2];
rz(-2.006532) q[2];
rz(0.48842397) q[3];
sx q[3];
rz(-1.6596183) q[3];
sx q[3];
rz(-1.9058913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9957073) q[0];
sx q[0];
rz(-1.1111525) q[0];
sx q[0];
rz(2.1671894) q[0];
rz(-2.3459332) q[1];
sx q[1];
rz(-2.7641422) q[1];
sx q[1];
rz(0.61765751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5159588) q[0];
sx q[0];
rz(-0.083481073) q[0];
sx q[0];
rz(0.030103695) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6154613) q[2];
sx q[2];
rz(-1.2210033) q[2];
sx q[2];
rz(-2.443231) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5646334) q[1];
sx q[1];
rz(-1.8398713) q[1];
sx q[1];
rz(-0.58522018) q[1];
rz(-pi) q[2];
rz(0.97891221) q[3];
sx q[3];
rz(-0.55044658) q[3];
sx q[3];
rz(1.5959671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.016102942) q[2];
sx q[2];
rz(-1.8718655) q[2];
sx q[2];
rz(1.7529091) q[2];
rz(-1.3056508) q[3];
sx q[3];
rz(-0.7235705) q[3];
sx q[3];
rz(1.2828264) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0131123) q[0];
sx q[0];
rz(-0.73902577) q[0];
sx q[0];
rz(2.541743) q[0];
rz(-2.5685617) q[1];
sx q[1];
rz(-2.532798) q[1];
sx q[1];
rz(2.1175687) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065627873) q[0];
sx q[0];
rz(-0.5813404) q[0];
sx q[0];
rz(-0.79728787) q[0];
rz(3.0863831) q[2];
sx q[2];
rz(-0.6533567) q[2];
sx q[2];
rz(-0.094815985) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3249859) q[1];
sx q[1];
rz(-1.4868951) q[1];
sx q[1];
rz(2.744426) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38174601) q[3];
sx q[3];
rz(-1.7110398) q[3];
sx q[3];
rz(2.563917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.92274252) q[2];
sx q[2];
rz(-1.7538193) q[2];
sx q[2];
rz(-0.77524033) q[2];
rz(-1.5357337) q[3];
sx q[3];
rz(-1.1865059) q[3];
sx q[3];
rz(1.2576013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6105462) q[0];
sx q[0];
rz(-2.1903867) q[0];
sx q[0];
rz(1.2280986) q[0];
rz(-1.3950521) q[1];
sx q[1];
rz(-1.4735305) q[1];
sx q[1];
rz(2.7458618) q[1];
rz(-0.45356815) q[2];
sx q[2];
rz(-1.8487052) q[2];
sx q[2];
rz(-2.1089274) q[2];
rz(1.6041605) q[3];
sx q[3];
rz(-2.3725474) q[3];
sx q[3];
rz(-0.65548246) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

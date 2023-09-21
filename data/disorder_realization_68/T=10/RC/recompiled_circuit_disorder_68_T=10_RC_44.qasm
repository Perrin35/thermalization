OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(2.6497901) q[0];
sx q[0];
rz(9.2368035) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(-1.6245276) q[1];
sx q[1];
rz(-2.7741073) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7029019) q[0];
sx q[0];
rz(-0.54649788) q[0];
sx q[0];
rz(1.370907) q[0];
x q[1];
rz(3.0066178) q[2];
sx q[2];
rz(-2.0581323) q[2];
sx q[2];
rz(-0.48263532) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5012813) q[1];
sx q[1];
rz(-0.2959364) q[1];
sx q[1];
rz(-2.179115) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3931307) q[3];
sx q[3];
rz(-1.9068309) q[3];
sx q[3];
rz(-0.48776585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1774896) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(-1.3163542) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47857639) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(2.6696894) q[0];
rz(0.42981237) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(0.93634161) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72915709) q[0];
sx q[0];
rz(-1.4834187) q[0];
sx q[0];
rz(2.9108414) q[0];
rz(-pi) q[1];
rz(-0.78511946) q[2];
sx q[2];
rz(-1.2650507) q[2];
sx q[2];
rz(2.608125) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97570005) q[1];
sx q[1];
rz(-2.4411538) q[1];
sx q[1];
rz(-0.16209929) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5283269) q[3];
sx q[3];
rz(-0.50308933) q[3];
sx q[3];
rz(2.344362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3669746) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(-2.7152087) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(-3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8957829) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(2.202503) q[0];
rz(-2.242873) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(-2.5476707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9868601) q[0];
sx q[0];
rz(-1.2504471) q[0];
sx q[0];
rz(2.8264168) q[0];
rz(3.0837768) q[2];
sx q[2];
rz(-2.2584492) q[2];
sx q[2];
rz(-0.85580326) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3214026) q[1];
sx q[1];
rz(-0.47649511) q[1];
sx q[1];
rz(-2.0501775) q[1];
rz(-2.6731554) q[3];
sx q[3];
rz(-1.2611946) q[3];
sx q[3];
rz(-2.7999511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5014191) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(-1.7017986) q[2];
rz(-0.38763186) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(-2.7534527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(-2.638812) q[0];
rz(2.373383) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(2.3847413) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605646) q[0];
sx q[0];
rz(-1.2111944) q[0];
sx q[0];
rz(-1.2147551) q[0];
x q[1];
rz(-3.1286131) q[2];
sx q[2];
rz(-1.105096) q[2];
sx q[2];
rz(-2.4109858) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4591603) q[1];
sx q[1];
rz(-1.9285413) q[1];
sx q[1];
rz(2.1898502) q[1];
rz(-pi) q[2];
rz(-0.54550708) q[3];
sx q[3];
rz(-2.3295998) q[3];
sx q[3];
rz(0.10520392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42671529) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(1.654401) q[2];
rz(2.5590844) q[3];
sx q[3];
rz(-1.094386) q[3];
sx q[3];
rz(2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-0.75772444) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(2.0910738) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3462853) q[0];
sx q[0];
rz(-1.2588358) q[0];
sx q[0];
rz(-2.8739268) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97425766) q[2];
sx q[2];
rz(-1.1968808) q[2];
sx q[2];
rz(2.8982382) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1041303) q[1];
sx q[1];
rz(-0.81172746) q[1];
sx q[1];
rz(2.3292259) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80661185) q[3];
sx q[3];
rz(-1.1685373) q[3];
sx q[3];
rz(-0.81392399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22333764) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(2.5081432) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(-0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(-0.25318405) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(-1.4621428) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40186858) q[0];
sx q[0];
rz(-1.7812294) q[0];
sx q[0];
rz(1.6519288) q[0];
rz(-pi) q[1];
rz(1.2573104) q[2];
sx q[2];
rz(-0.67337155) q[2];
sx q[2];
rz(-1.9999258) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.05789214) q[1];
sx q[1];
rz(-1.18827) q[1];
sx q[1];
rz(1.7544569) q[1];
x q[2];
rz(1.2268279) q[3];
sx q[3];
rz(-1.2729537) q[3];
sx q[3];
rz(2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2016466) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(-2.3366826) q[2];
rz(-1.9645875) q[3];
sx q[3];
rz(-1.0369119) q[3];
sx q[3];
rz(2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2729623) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(1.0246798) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(2.129508) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9818078) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(0.83321379) q[0];
rz(-pi) q[1];
rz(-0.32832844) q[2];
sx q[2];
rz(-0.40847455) q[2];
sx q[2];
rz(2.7834053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.583657) q[1];
sx q[1];
rz(-0.54469889) q[1];
sx q[1];
rz(-3.0779561) q[1];
rz(-pi) q[2];
rz(0.55862553) q[3];
sx q[3];
rz(-1.6405676) q[3];
sx q[3];
rz(-2.74673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6107789) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(0.42993316) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(-0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5291418) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(0.21417831) q[0];
rz(-1.051349) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(-2.8578551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9203628) q[0];
sx q[0];
rz(-1.2698679) q[0];
sx q[0];
rz(-1.5350585) q[0];
rz(-pi) q[1];
rz(1.9010504) q[2];
sx q[2];
rz(-1.7527765) q[2];
sx q[2];
rz(0.44588003) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2092065) q[1];
sx q[1];
rz(-1.3750409) q[1];
sx q[1];
rz(0.50435658) q[1];
rz(-pi) q[2];
rz(0.78332087) q[3];
sx q[3];
rz(-2.3187227) q[3];
sx q[3];
rz(0.93572039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45067898) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(1.4452176) q[2];
rz(-1.5971659) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-2.8022695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63672367) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(3.0723363) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5690631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78281392) q[0];
sx q[0];
rz(-0.94952119) q[0];
sx q[0];
rz(-0.57399477) q[0];
rz(-pi) q[1];
rz(1.6530767) q[2];
sx q[2];
rz(-1.2459718) q[2];
sx q[2];
rz(-2.5282853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.931817) q[1];
sx q[1];
rz(-1.6546384) q[1];
sx q[1];
rz(-0.2340338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2194355) q[3];
sx q[3];
rz(-1.4033917) q[3];
sx q[3];
rz(0.64099121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1853603) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(3.0017079) q[2];
rz(2.774003) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(-0.99115133) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-0.26783255) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5726686) q[0];
sx q[0];
rz(-2.2974282) q[0];
sx q[0];
rz(2.5323244) q[0];
rz(-3.015976) q[2];
sx q[2];
rz(-1.3889424) q[2];
sx q[2];
rz(2.6806718) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2512868) q[1];
sx q[1];
rz(-2.7000513) q[1];
sx q[1];
rz(2.4114386) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0288826) q[3];
sx q[3];
rz(-1.3241323) q[3];
sx q[3];
rz(2.820462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3315167) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(2.5951071) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(0.77990445) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(-0.71074625) q[2];
sx q[2];
rz(-2.032861) q[2];
sx q[2];
rz(-2.1326333) q[2];
rz(0.6775425) q[3];
sx q[3];
rz(-1.4031706) q[3];
sx q[3];
rz(-1.8085898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
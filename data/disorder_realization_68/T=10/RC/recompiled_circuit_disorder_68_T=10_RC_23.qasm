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
rz(-0.49180254) q[0];
sx q[0];
rz(0.1879745) q[0];
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(-0.36748537) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039283218) q[0];
sx q[0];
rz(-1.4674205) q[0];
sx q[0];
rz(1.0331868) q[0];
rz(-pi) q[1];
rz(-1.3221402) q[2];
sx q[2];
rz(-0.50422943) q[2];
sx q[2];
rz(0.76489514) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6403113) q[1];
sx q[1];
rz(-0.2959364) q[1];
sx q[1];
rz(2.179115) q[1];
x q[2];
rz(1.3931307) q[3];
sx q[3];
rz(-1.2347617) q[3];
sx q[3];
rz(-2.6538268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.964103) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(-2.5906079) q[2];
rz(-1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.47857639) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(-0.4719032) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(-0.93634161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72915709) q[0];
sx q[0];
rz(-1.658174) q[0];
sx q[0];
rz(-0.23075128) q[0];
rz(-0.41994862) q[2];
sx q[2];
rz(-0.83050767) q[2];
sx q[2];
rz(-2.3967957) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1658926) q[1];
sx q[1];
rz(-2.4411538) q[1];
sx q[1];
rz(-2.9794934) q[1];
x q[2];
rz(-1.5283269) q[3];
sx q[3];
rz(-0.50308933) q[3];
sx q[3];
rz(-2.344362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77461809) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(-1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(-0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24580978) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(-2.202503) q[0];
rz(2.242873) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(2.5476707) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.827841) q[0];
sx q[0];
rz(-1.272164) q[0];
sx q[0];
rz(1.9065501) q[0];
rz(-pi) q[1];
rz(0.8823231) q[2];
sx q[2];
rz(-1.6154628) q[2];
sx q[2];
rz(2.4633212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8507104) q[1];
sx q[1];
rz(-1.1516654) q[1];
sx q[1];
rz(0.23371975) q[1];
rz(-pi) q[2];
rz(1.2265615) q[3];
sx q[3];
rz(-1.1262745) q[3];
sx q[3];
rz(1.759474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5014191) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(-1.4397941) q[2];
rz(-0.38763186) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(2.7534527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7664292) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(2.638812) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(-2.3847413) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2198974) q[0];
sx q[0];
rz(-1.9031525) q[0];
sx q[0];
rz(0.38145782) q[0];
rz(-0.012979522) q[2];
sx q[2];
rz(-1.105096) q[2];
sx q[2];
rz(-0.73060689) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0087352) q[1];
sx q[1];
rz(-0.99616226) q[1];
sx q[1];
rz(2.7112714) q[1];
rz(-0.73369153) q[3];
sx q[3];
rz(-1.1847704) q[3];
sx q[3];
rz(-1.2804077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42671529) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(-1.654401) q[2];
rz(-2.5590844) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29397598) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-0.75772444) q[0];
rz(1.2879397) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(1.0505189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3462853) q[0];
sx q[0];
rz(-1.8827569) q[0];
sx q[0];
rz(-2.8739268) q[0];
x q[1];
rz(-0.97425766) q[2];
sx q[2];
rz(-1.9447118) q[2];
sx q[2];
rz(0.24335441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.95549059) q[1];
sx q[1];
rz(-2.0932066) q[1];
sx q[1];
rz(-2.2239457) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0195144) q[3];
sx q[3];
rz(-2.2973804) q[3];
sx q[3];
rz(2.7725078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.22333764) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(-0.63344947) q[2];
rz(-1.194362) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(2.4244394) q[3];
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
rz(-0.69960064) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(2.8884086) q[0];
rz(1.6075915) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(1.4621428) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40186858) q[0];
sx q[0];
rz(-1.3603633) q[0];
sx q[0];
rz(-1.4896638) q[0];
rz(2.2199549) q[2];
sx q[2];
rz(-1.3772794) q[2];
sx q[2];
rz(0.18093872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51965442) q[1];
sx q[1];
rz(-2.7192273) q[1];
sx q[1];
rz(-2.7154891) q[1];
rz(1.2268279) q[3];
sx q[3];
rz(-1.868639) q[3];
sx q[3];
rz(-2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2016466) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(2.3366826) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2729623) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(-0.92765635) q[0];
rz(-1.0246798) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(-1.0120846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0976809) q[0];
sx q[0];
rz(-1.0543531) q[0];
sx q[0];
rz(0.54215706) q[0];
rz(-pi) q[1];
rz(-1.7094678) q[2];
sx q[2];
rz(-1.1853293) q[2];
sx q[2];
rz(-3.1388381) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.583657) q[1];
sx q[1];
rz(-0.54469889) q[1];
sx q[1];
rz(3.0779561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6530232) q[3];
sx q[3];
rz(-2.1279018) q[3];
sx q[3];
rz(1.132387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6107789) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(2.7116595) q[2];
rz(-2.1271465) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(-0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5291418) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(-2.9274143) q[0];
rz(1.051349) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(0.28373757) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9203628) q[0];
sx q[0];
rz(-1.2698679) q[0];
sx q[0];
rz(1.6065341) q[0];
rz(1.2405422) q[2];
sx q[2];
rz(-1.3888161) q[2];
sx q[2];
rz(0.44588003) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2092065) q[1];
sx q[1];
rz(-1.7665518) q[1];
sx q[1];
rz(-0.50435658) q[1];
x q[2];
rz(-2.22105) q[3];
sx q[3];
rz(-2.1170108) q[3];
sx q[3];
rz(-1.2342681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45067898) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(-1.4452176) q[2];
rz(1.5444267) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504869) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(0.069256393) q[0];
rz(-1.6537369) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5725296) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3587787) q[0];
sx q[0];
rz(-0.94952119) q[0];
sx q[0];
rz(-0.57399477) q[0];
rz(0.23937978) q[2];
sx q[2];
rz(-0.33472543) q[2];
sx q[2];
rz(2.7811546) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3410586) q[1];
sx q[1];
rz(-1.8039928) q[1];
sx q[1];
rz(1.6569767) q[1];
x q[2];
rz(1.2980372) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(2.42815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1853603) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(0.13988477) q[2];
rz(2.774003) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(-0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-0.26783255) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7627015) q[0];
sx q[0];
rz(-0.91103103) q[0];
sx q[0];
rz(0.99878175) q[0];
rz(-pi) q[1];
rz(0.97271131) q[2];
sx q[2];
rz(-0.22062606) q[2];
sx q[2];
rz(2.0711183) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2512868) q[1];
sx q[1];
rz(-0.44154134) q[1];
sx q[1];
rz(-2.4114386) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0288826) q[3];
sx q[3];
rz(-1.8174603) q[3];
sx q[3];
rz(-0.32113069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(1.4769185) q[2];
rz(2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(-0.5464856) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289566) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(0.77990445) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(-2.1521679) q[2];
sx q[2];
rz(-0.94716723) q[2];
sx q[2];
rz(-0.92826044) q[2];
rz(-0.26364636) q[3];
sx q[3];
rz(-2.4468138) q[3];
sx q[3];
rz(3.1082603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

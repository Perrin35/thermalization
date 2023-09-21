OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(0.89755091) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7449887) q[0];
sx q[0];
rz(-0.94192266) q[0];
sx q[0];
rz(0.086358503) q[0];
rz(-0.35383309) q[2];
sx q[2];
rz(-2.1733123) q[2];
sx q[2];
rz(1.8941855) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7251016) q[1];
sx q[1];
rz(-0.93826586) q[1];
sx q[1];
rz(-1.3593258) q[1];
rz(-pi) q[2];
rz(-0.75818054) q[3];
sx q[3];
rz(-1.4361793) q[3];
sx q[3];
rz(-1.4338223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66427461) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(2.409639) q[2];
rz(-0.96015635) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(-1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(0.91645855) q[0];
rz(-0.48049277) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(-2.2629471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4449578) q[0];
sx q[0];
rz(-2.1793471) q[0];
sx q[0];
rz(0.47136013) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36388134) q[2];
sx q[2];
rz(-0.65081396) q[2];
sx q[2];
rz(-1.3712856) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6138184) q[1];
sx q[1];
rz(-1.6442181) q[1];
sx q[1];
rz(-0.41927494) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7813803) q[3];
sx q[3];
rz(-0.31616351) q[3];
sx q[3];
rz(-2.8663243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(-0.64727616) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6140401) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(0.24599427) q[0];
rz(1.4100769) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(1.1211959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778053) q[0];
sx q[0];
rz(-2.2664547) q[0];
sx q[0];
rz(1.0870766) q[0];
rz(0.063935117) q[2];
sx q[2];
rz(-1.1786412) q[2];
sx q[2];
rz(-2.8901697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9417291) q[1];
sx q[1];
rz(-0.45048303) q[1];
sx q[1];
rz(-0.73168879) q[1];
x q[2];
rz(-1.6882012) q[3];
sx q[3];
rz(-2.061764) q[3];
sx q[3];
rz(1.9672058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-2.1901954) q[2];
rz(0.65008632) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-0.78805584) q[0];
rz(-1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1040092) q[0];
sx q[0];
rz(-1.5452256) q[0];
sx q[0];
rz(2.5081162) q[0];
x q[1];
rz(-2.7599081) q[2];
sx q[2];
rz(-1.0978062) q[2];
sx q[2];
rz(2.4406529) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1981922) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(2.0342159) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4898407) q[3];
sx q[3];
rz(-1.6162989) q[3];
sx q[3];
rz(-0.10304606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5300166) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(-2.7799515) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5761121) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(-0.53260032) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(-1.9929569) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38265739) q[0];
sx q[0];
rz(-1.4969345) q[0];
sx q[0];
rz(-2.1942684) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0257452) q[2];
sx q[2];
rz(-0.75136853) q[2];
sx q[2];
rz(0.92912208) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7733113) q[1];
sx q[1];
rz(-1.9241663) q[1];
sx q[1];
rz(1.7721304) q[1];
x q[2];
rz(-1.9653737) q[3];
sx q[3];
rz(-2.2567281) q[3];
sx q[3];
rz(1.3057614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0107515) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(0.71470913) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(0.39598879) q[0];
rz(-1.4453325) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(-0.2063624) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43564046) q[0];
sx q[0];
rz(-2.3898976) q[0];
sx q[0];
rz(-3.0138426) q[0];
rz(2.2588737) q[2];
sx q[2];
rz(-2.1157017) q[2];
sx q[2];
rz(-2.8628652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3316162) q[1];
sx q[1];
rz(-0.14252256) q[1];
sx q[1];
rz(-0.82377164) q[1];
x q[2];
rz(1.2431074) q[3];
sx q[3];
rz(-1.7565691) q[3];
sx q[3];
rz(-2.1873308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.53987327) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(-2.8708141) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8966184) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(-1.2111838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9528708) q[0];
sx q[0];
rz(-2.0824008) q[0];
sx q[0];
rz(0.94922025) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6518875) q[2];
sx q[2];
rz(-1.4486794) q[2];
sx q[2];
rz(-0.44039886) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3019575) q[1];
sx q[1];
rz(-2.9710899) q[1];
sx q[1];
rz(-1.9819928) q[1];
rz(-0.32254036) q[3];
sx q[3];
rz(-0.99223677) q[3];
sx q[3];
rz(1.6998147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3367735) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(-0.024519196) q[2];
rz(-0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(1.0205644) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.9205836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97604254) q[0];
sx q[0];
rz(-1.3001469) q[0];
sx q[0];
rz(-0.793215) q[0];
x q[1];
rz(-1.4079354) q[2];
sx q[2];
rz(-1.3961627) q[2];
sx q[2];
rz(2.3101431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11215969) q[1];
sx q[1];
rz(-2.9661848) q[1];
sx q[1];
rz(2.4584241) q[1];
rz(1.9302759) q[3];
sx q[3];
rz(-0.7822789) q[3];
sx q[3];
rz(-0.6560916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(-2.4364046) q[2];
rz(-2.5382036) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(2.5367472) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(0.12577122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98709479) q[0];
sx q[0];
rz(-0.25941601) q[0];
sx q[0];
rz(0.82018606) q[0];
x q[1];
rz(0.33340402) q[2];
sx q[2];
rz(-1.7575022) q[2];
sx q[2];
rz(-2.8508027) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3667664) q[1];
sx q[1];
rz(-0.80087304) q[1];
sx q[1];
rz(-1.1714539) q[1];
rz(-pi) q[2];
rz(1.3689234) q[3];
sx q[3];
rz(-1.2807506) q[3];
sx q[3];
rz(1.2253075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13828364) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9897292) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.9412769) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(-2.0013924) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32894293) q[0];
sx q[0];
rz(-2.86033) q[0];
sx q[0];
rz(-1.2055231) q[0];
rz(-pi) q[1];
rz(-2.0613725) q[2];
sx q[2];
rz(-0.51418257) q[2];
sx q[2];
rz(1.0487923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0298437) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(-2.7010121) q[1];
rz(1.2530008) q[3];
sx q[3];
rz(-2.3386049) q[3];
sx q[3];
rz(-0.12249891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
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
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-0.76277914) q[2];
sx q[2];
rz(-0.63763466) q[2];
sx q[2];
rz(0.55479738) q[2];
rz(1.981001) q[3];
sx q[3];
rz(-2.1174178) q[3];
sx q[3];
rz(-2.5934364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
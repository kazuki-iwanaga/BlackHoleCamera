FROM python:3.7.5-alpine3.10

RUN apk update && \
    apk add --no-cache gcc g++ gfortran \
    libc-dev boost-dev musl-dev lapack lapack-dev make cmake git

WORKDIR /root/myenv

CMD [ "/bin/sh" ]
